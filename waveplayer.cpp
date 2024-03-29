#include <iostream>
#include <cmath>
#include <Windows.h>
#include <portaudio.h>
#include "misc.h"

#define SAMPLE_RATE   (44100)
#define FRAMES_PER_BUFFER  (256)
#define FREQUENCY (440)

#ifndef M_PI
#define M_PI  (3.14159265)
#endif

#define TABLE_SIZE   (200)

using namespace std;

struct Envelope
{
	double dAttackTime;
	double dDecayTime;
	double dSustainAmplitude;
	double dReleaseTime;
	double dStartAmplitude;
	double dTriggerOffTime;
	double dTriggerOnTime;
	bool bNoteOn;

	Envelope()
	{
		dAttackTime = 0.10;
		dDecayTime = 0.02;
		dStartAmplitude = 1.0;
		dSustainAmplitude = 0.8;
		dReleaseTime = 2.0;
		bNoteOn = false;
		dTriggerOffTime = 0.0;
		dTriggerOnTime = 0.0;
	}

	void NoteOn(double dTimeOn)
	{
		dTriggerOnTime = dTimeOn;
		bNoteOn = true;
	}

	void NoteOff(double dTimeOff)
	{
		dTriggerOffTime = dTimeOff;
		bNoteOn = false;
	}

	double GetAmplitude(double dTime)
	{
		double dAmplitude = 0.0;
		double dLifeTime = dTime - dTriggerOnTime;

		if (bNoteOn)
		{
			if (dLifeTime <= dAttackTime)
			{
				// In attack Phase - approach max amplitude
				dAmplitude = (dLifeTime / dAttackTime) * dStartAmplitude;
			}

			if (dLifeTime > dAttackTime && dLifeTime <= (dAttackTime + dDecayTime))
			{
				// In decay phase - reduce to sustained amplitude
				dAmplitude = ((dLifeTime - dAttackTime) / dDecayTime) * (dSustainAmplitude - dStartAmplitude) + dStartAmplitude;
			}

			if (dLifeTime > (dAttackTime + dDecayTime))
			{
				// In sustain phase - dont change until note released
				dAmplitude = dSustainAmplitude;
			}
		}
		else
		{
			// Note has been released, so in release phase
			dAmplitude = ((dTime - dTriggerOffTime) / dReleaseTime) * (0.0 - dSustainAmplitude) + dSustainAmplitude;
		}

		// Amplitude should not be negative
		if (dAmplitude <= 0.0001)
			dAmplitude = 0.0;

		return dAmplitude;
	}
};


int main()
{
	PaStreamParameters outputParameters;
	PaStream *stream;
	PaError err;
	float buffer[FRAMES_PER_BUFFER][2];
	GTABLE* gtable;
	OSCILT* oscilt;
	int i, j;
	int buffercount;
	double freq;
	double basefreq;
	double d12thRootOf2;
	Envelope env;
	double amp;
	bool keyPressed;
	double time, timeincr;
	int currentkey;

	gtable = new_saw(TABLE_SIZE, 5, 0);
	oscilt = new_oscilt(SAMPLE_RATE, gtable, 0);

	err = Pa_Initialize();
	if (err != paNoError)
		goto error;

	outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
	if (outputParameters.device == paNoDevice)
	{
		cout << "Error: No default output device.\n";
		goto error;
	}
	outputParameters.channelCount = 2;       /* stereo output */
	outputParameters.sampleFormat = paFloat32; /* 32 bit floating point output */
	outputParameters.suggestedLatency = Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
	outputParameters.hostApiSpecificStreamInfo = NULL;

	err = Pa_OpenStream(
		&stream,
		NULL, /* no input */
		&outputParameters,
		SAMPLE_RATE,
		FRAMES_PER_BUFFER,
		paClipOff,      /* we won't output out of range samples so don't bother clipping them */
		NULL,
		NULL);
	if (err != paNoError)
		goto error;

	err = Pa_StartStream(stream);
	if (err != paNoError)
		goto error;

	buffercount = 8;
	basefreq = 440;
	freq = 0;
	d12thRootOf2 = pow(2.0, 1.0 / 12.0);
	timeincr = 1.0 / SAMPLE_RATE;
	time = 0;
	currentkey = -1;

	while (1)
	{
		keyPressed = false;
		for (int k = 0; k < 16; k++)
		{
			if (GetAsyncKeyState((unsigned char)("ZSXCFVGBNJMK\xbcL\xbe\xbf"[k])) & 0x8000)
			{
				if (currentkey != k)
				{
					freq = basefreq * pow(d12thRootOf2, k);
					env.NoteOn(time);
					cout << "\rNote On : " << time << "s " << freq << "Hz";
					currentkey = k;
				}
				keyPressed = true;
			}
		}
		if(!keyPressed)
		{
			if (currentkey != -1)
			{
				cout << "\rNote Off: " << time << "s";
				env.NoteOff(time);
				currentkey = -1;
			}
		}
		
		for (i = 0; i < buffercount; i++)
		{
			for (j = 0; j < FRAMES_PER_BUFFER; j++)
			{
				buffer[j][0] = env.GetAmplitude(time) * float(tick(oscilt, freq));
				buffer[j][1] = buffer[j][0];
				time += timeincr;
			}
			err = Pa_WriteStream(stream, buffer, FRAMES_PER_BUFFER);
			if (err != paNoError)
				goto error;
		}
		if (GetAsyncKeyState(VK_UP))
			break;
	}
	err = Pa_StopStream(stream);
	if (err != paNoError)
		goto error;

	err = Pa_CloseStream(stream);
	if (err != paNoError) goto error;

	Pa_Terminate();
	printf("Test finished.\n");

	return err;

error:
	fprintf(stderr, "An error occured while using the portaudio stream\n");
	fprintf(stderr, "Error number: %d\n", err);
	fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(err));
	// Print more information about the error.
	if (err == paUnanticipatedHostError)
	{
		const PaHostErrorInfo *hostErrorInfo = Pa_GetLastHostErrorInfo();
		fprintf(stderr, "Host API error = #%ld, hostApiType = %d\n", hostErrorInfo->errorCode, hostErrorInfo->hostApiType);
		fprintf(stderr, "Host API error = %s\n", hostErrorInfo->errorText);
	}
	Pa_Terminate();
	return err;
}
