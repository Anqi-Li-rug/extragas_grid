#ifndef FINDVALUE_H
#define FINDVALUE_H

#include <iostream>
#include <cstdlib>
#include "findValue.h"
#include <stdio.h>
#include <string.h>
void findValue(char filename[], int &value, char parameter[])
{
  FILE* pinput_file;
  pinput_file=fopen(filename,"r");
  char par_name[80], par_value[80];
  while (strncmp(par_name, parameter, strlen(parameter))) {
    fscanf (pinput_file, "%s ", &par_name);
    fscanf (pinput_file, "%s ", &par_value);
  }
  value=atoi(par_value);
  fclose(pinput_file);
}

void findValue(char filename[], float &value, char parameter[])
{
  FILE* pinput_file;
  pinput_file=fopen(filename,"r");
  char par_name[80], par_value[80];
  while (strncmp(par_name, parameter, strlen(parameter))) {
    fscanf (pinput_file, "%s ", &par_name);
    fscanf (pinput_file, "%s ", &par_value);
  }
  value=atof(par_value);
  fclose(pinput_file);
}

void findValue(char filename[], double &value, char parameter[])
{
  FILE* pinput_file;
  pinput_file=fopen(filename,"r");
  char par_name[80], par_value[80];
  while (strncmp(par_name, parameter, strlen(parameter))) {
    fscanf (pinput_file, "%s ", &par_name);
    fscanf (pinput_file, "%s ", &par_value);
  }
  value=(double)atof(par_value);
  fclose(pinput_file);
}

void findValue(char filename[], char value[], char parameter[])
{
  FILE* pinput_file;
  pinput_file=fopen(filename,"r");
  char par_name[80], par_value[80];
  while (strncmp(par_name, parameter, strlen(parameter))) {
    fscanf (pinput_file, "%s ", &par_name);
    fscanf (pinput_file, "%s ", &par_value);
  }
  strcpy(value,par_value);
  fclose(pinput_file);
}

void findValue(char filename[], char &value, char parameter[])
{
  FILE* pinput_file;
  pinput_file=fopen(filename,"r");
  char par_name[80], par_value[80];
  while (strncmp(par_name, parameter, strlen(parameter))) {
    fscanf (pinput_file, "%s ", &par_name);
    fscanf (pinput_file, "%s ", &par_value);
  }
  value=par_value[0];
  fclose(pinput_file);
}

void findValue(char filename[], bool &value, char parameter[])
{
	FILE* pinput_file;
	pinput_file=fopen(filename,"r");
	char par_name[80], par_value[80];
	while (strncmp(par_name, parameter, strlen(parameter))) {
		fscanf (pinput_file, "%s ", &par_name);
		fscanf (pinput_file, "%s ", &par_value);
	}
	if(par_value[0]=='y') {value=true;}
	else {value=false;}
	fclose(pinput_file);
}
#endif
