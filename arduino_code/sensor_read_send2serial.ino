#include "DFRobot_AS7341.h"
#include <Wire.h>        
#include <Arduino.h>
#include <SPI.h>
#include <stdio.h>
#include <stdlib.h>

DFRobot_AS7341 as7341;

byte adr = 0x39; // Sensor address

float F1 = 0, F2 = 0, F3 = 0, F4 = 0, F5 = 0, F6 = 0, F7 = 0, F8 = 0, NIR = 0, clear = 0, flicker = 0;

const uint8_t ATIME = 0;
const uint16_t ASTEP = 65534;
const uint16_t GAIN512 = 512; // For F1 - F6
const uint16_t GAIN256 = 256; // For F7, F8, Clear

float tint_ms = 0;

void setup() {
  Serial.begin(9600);
  as7341.begin();

  /* Time configurations -------------------------------------------------------------------------------------------------*/
  //Set the value of register ATIME, through which the value of Integration time can be calculated. The value represents the time that must be spent during data reading.
  as7341.setAtime(ATIME);

  //Set the value of register ASTEP, through which the value of Integration time can be calculated. The value represents the time that must be spent during data reading.
  as7341.setAstep(ASTEP);

  //Integration time = (ATIME + 1) x (ASTEP + 1) x 2.78µs
  tint_ms = (ATIME + 1) * (ASTEP + 1) * (2.78/1000);
}

void loop() {
// Run sensor code

  DFRobot_AS7341::sModeOneData_t data1;
  DFRobot_AS7341::sModeTwoData_t data2;

  as7341.setAGAIN(10);

  //Start spectrum measurement 
  //Channel mapping mode: 1.eF1F4ClearNIR,2.eF5F8ClearNIR
  as7341.startMeasure(as7341.eF1F4ClearNIR);
  
  //Read the value of sensor data channel 0~5, under eF1F4ClearNIR
  data1 = as7341.readSpectralDataOne();

  // BasicCounts = (RawCounts - OffSet) / (gain * tint_ms)
  //F1 = (data1.ADF1 / (GAIN512 * tint_ms));
  Serial.println(data1.ADF1);

  //F2 = (data1.ADF2 / (GAIN512 * tint_ms));
  Serial.println(data1.ADF2);

  //F3 = (data1.ADF3 / (GAIN512 * tint_ms));
  Serial.println(data1.ADF3);

  //F4 = (data1.ADF4 / (GAIN512 * tint_ms));
  Serial.println(data1.ADF4);
  
  as7341.startMeasure(as7341.eF5F8ClearNIR);    
  data2 = as7341.readSpectralDataTwo();    

  //F5 = (data2.ADF5 / (GAIN512 * tint_ms));
  Serial.println(data2.ADF5);

  //F6 = (data2.ADF6 / (GAIN512 * tint_ms));
  Serial.println(data2.ADF6);

  //Set gain value(0~10 corresponds to X0.5,X1,X2,X4,X8,X16,X32,X64,X128,X256,X512)
  as7341.setAGAIN(9);

  // Get sensor values again because the gain for F7, F8, and clear are different
  as7341.startMeasure(as7341.eF5F8ClearNIR);    
  data2 = as7341.readSpectralDataTwo(); 

  //F7 = (data2.ADF7 / (GAIN256 * tint_ms));
  Serial.println(data2.ADF7);

  //F8 = (data2.ADF8 / (GAIN256 * tint_ms));
  Serial.println(data2.ADF8);

  //clear = (data2.ADCLEAR / (GAIN256 * tint_ms));
  Serial.println(data2.ADCLEAR);
  
  //NIR = (data2.ADNIR / (GAIN256 * tint_ms));
  Serial.println(data2.ADNIR);
  
  Wire.endTransmission();
  delay(10);
}