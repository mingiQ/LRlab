/*  Leonid Rokhinson © 2022  MicroLambda BPF type ### control  
*/
#include <SPI.h>
#include <stdio.h>
#include <strings.h>
#include <EEPROM.h>

#define devstamp  "uLambda BPF v1.0 (c) 11-25-2022"

#define maxmem 360000             // max dinamic memory ("... leaving XXXXXXX bytes for local variables")

/*    ------ communication commands ---------
 *     serial @ 115200, line end character: \n
 *    "?" or "*IDN?"             ->  devstamp         - equipment ID  
 *    "FF val"                   -> "FF"              - set single dac value (float)
 *    "SF fmin fmax"             -> "SF"              - save fmin and fmax values to EEPROM
 *    "RF"                       -> "RF"              - read fmin and fmax values from EEPROM
 *    
 *    "TS dac dval"              -> "TS"              - comm test: write xAAAAA to dac dval times
 *    "TX pin dval"              -> "TX"              - pin test: write HIGH/LOW to pin dval times
 *    "TT pin"                   -> "TT"              - toggle pin
*/

#define CS           10       // chip select

#define spi_MOSI     11       // dac_SDO
#define spi_MISO     12       // dac_SDIN
#define spi_SCLK     13       // data on high-to-low, idle high

#define spi_speed   5000000    // speed in Hz 5MHz
#define spi_bit     MSBFIRST    // 
#define spi_mode    SPI_MODE0   // SCLK idle LOW, data capture low-to-high
unsigned char spibuff[2];       // buffer for SPI write
//#define spi_mode    SPI_MODE3   // on T3.6? SCLK idle HIGH, data capture low-to-high
// SPI.beginTransaction() & SPI.endTransaction() do we need to release the SPI ?


#define ee_start_addr 0x00    // EEPROM address to store dacpars
#define pi2     6.28318530718 // 2pi
#define lb20       0x7FFFF    // 2^19-1
#define lb18       0x1FFFF    // 2^17-1

struct EEPROM_t {             // EEPROM record
  float fqmin;
  float fqmax;
};

String  inputString = "";                   // a string to hold command
volatile boolean stringComplete = false;    // whether the string is complete
volatile boolean IRQready;                  // flag processed via serial interrupt  
volatile boolean cmp_flag;                  // flag processed by comparator interrupt

String inputBuffer = "";                    // a string to hold incoming data
char buff[256];                             // string buffer 
char scommand[64];                          // string to hold command (at least TFT-wide)
uint32_t comtimer;                          // counts ms for the input transmission


unsigned int ibts, lval, pin;
float freq,fqmin,fqmax;

elapsedMicros wtime, ustimer;               // us timer
uint32_t loop_time;                         // us for the last loop
uint32_t t,tt,tmax,tmin=10000;              // vars for us counter and operations
String stamp;                               // version stamp
boolean update_flag;                        // whether update LDAC, ones per loop

EEPROM_t ee_saved, ee_new;                  // record to write to EEPROM
uint16_t ee_addr;

elapsedMicros usec;                            // for IC,IR commands

void setup() {
//  mbuff  = (float*) malloc ( maxmem );         // allocate dinamic memory for temp storage
//  mmbuff = (float*) malloc ( 4 );              // address of the last mbuff for debugs
//  maxmemfloats = maxmem / sizeof(float);       
  
  pinMode(CS, OUTPUT);
 
  digitalWrite(CS, HIGH);            // CS idle HIGH
  Serial.begin(115200);  
  Serial.setTimeout(1);              // initialize serial

  inputString.reserve(256);          // reserve bytes for the inputString
  stamp=String(String(devstamp)); //  + String(tver));

  SPI.begin();                            // initialize SPI for dacs

//  EEPROM.get(ee_start_addr, ee_saved);
//  fqmax = ee_saved.fqmax;
//  fqmin = ee_saved.fqmin;
  fqmax = 2000000000;
  fqmin = 18000000000;

  IRQready = false;
  cmp_flag = false;
}

void loop() {
 // Serial.println(" Hello from teensy");
//  update_flag = false;  // use COM port to communicate with arduino 2022-12-27
  while (!Serial.available());
  inputString = Serial.readString();
  inputString.toUpperCase();
  inputString.toCharArray(buff, inputString.length()+1);
  inputString = "";
  sscanf(buff, "%s", scommand);
  //if (stringComplete) {
  //  stringComplete=false;
   // inputString.toUpperCase();
   // inputString.toCharArray(buff, inputString.length()+1);
   // inputString = "";  
   // sscanf(buff, "%s", scommand);                                                     // read command from buff

/* ? & *IND? --- device ID responce  */
    if ( (strcmp(scommand, "?") == 0) || (strcmp(scommand, "*IDN?") == 0) ) {    
      Serial.println(stamp);
/* FF -- set frequency */
    } else if ( strcmp(scommand, "FF") == 0 ) {                          
      sscanf(&buff[2], "%g", &freq);
      lval = (freq-fqmin)/(fqmax-fqmin)*65535;
      write_SPI(lval);      // set freq
      sprintf(buff, "FF %d", lval);
      Serial.println(buff);   
/* FH -- set frequency */
    } else if ( strcmp(scommand, "FH") == 0 ) {                          
      sscanf(&buff[2], "%X", &lval);
      write_SPI(lval);      // set freq
      sprintf(buff, "FH 0x%X", lval);
      Serial.println(buff);   
/* SF -- write data to EEPROM --- */
    } else if ( strcmp(scommand, "SF") == 0 ) {                               // rewrites to EEPROM only if changed
        sscanf(&buff[2], "%g %g", &fqmin, &fqmax);
        ee_saved.fqmax = fqmax;
        ee_saved.fqmin = fqmin;
        EEPROM.put(ee_start_addr, ee_saved);  
        Serial.println("SF");
/* RF -- read data from EEPROM --- */
    } else if ( strcmp(scommand, "RF") == 0 ) {                               // rewrites to EEPROM only if changed
        EEPROM.get(ee_start_addr, ee_saved);
        fqmax = ee_saved.fqmax;
        fqmin = ee_saved.fqmin;
        sprintf(buff, "RF %f %f", fqmin, fqmax);
        Serial.println(buff);
/* TS -- multiple writes for DIO test --- */
    } else if ( strcmp(scommand, "TS") == 0 ) {                                 
      sscanf(&buff[2], "%u %x", &ibts, &lval);
        //lval = 0x0000;
        for(unsigned int i=0;i<ibts;i++){
          write_SPI(lval);
          delayMicroseconds(20);
        }
        Serial.println("TS");
/* TT -- toggle pin --- */
    } else if ( strcmp(scommand, "TT") == 0 ) {                                 
      sscanf(&buff[2], "%d", &pin);
      pinMode(pin, OUTPUT);
//      if (  digitalRead(ndac)==HIGH ) digitalWrite(pin, LOW );                        
//        else                          digitalWrite(pin, HIGH);                       
      Serial.println("TT");
/* TX -- writes to pin test --- */
    } else if ( strcmp(scommand, "TX") == 0 ) {                                 
      sscanf(&buff[2], "%d %u", &pin, &ibts);
      pinMode(pin, OUTPUT);
//      t=micros();
      for(unsigned int i=0;i<ibts;i++){
        digitalWrite(pin, LOW );                        
        digitalWrite(pin, HIGH);                       
      }
//      t=micros()-t;
      digitalWrite(pin, LOW );                        
      Serial.println("TX");
  }  else {
      Serial.println("2");   // wrong command
  }; /* end if ( scommand ) */
 //}; /* end if ( stringComplete )  */       
//  noInterrupts();  
//  t=micros();                                                                         // absolute us timer, resets every 11 min !!!
//  loop_time = wtime;                                                                  // us for the previous loop
//  wtime = 0;                                                                          // reset loop timer

}   // <<<<------------------------------ end loop() ------------------------------------->>>>>>>


/* write to SPI */
void write_SPI(unsigned int llval) {
  digitalWrite(CS, LOW);                             // enable input register
  spibuff[0] = (byte)((llval >> 8 ) & 0XFF);
  spibuff[1] = (byte)( llval        & 0XFF);
  SPI.beginTransaction(SPISettings(spi_speed, spi_bit, spi_mode));  // ?
  SPI.transfer(spibuff,2);
  SPI.endTransaction();                                             // ?
  digitalWrite(CS, HIGH);                            // update input register
}

/*  writes a HEX number to a single DAC AD57xx  
uint32_t  write_AD57xx_HEX(int dac, int32_t llval) {
  daci=dactable[dac]-1;
  llval = ( llval >  lmax ) ?  lmax : llval;
  llval = ( llval < -lmax ) ? -lmax : llval;        
  digitalWrite(board_addr_l, (daci >> 2) & 0x1 );                // select board
  digitalWrite(board_addr_m, (daci >> 3) & 0x1 );                // ----- || ----- 
  digitalWrite(board_addr_h, (daci >> 4) & 0x1 );                // ----- || -----
  digitalWrite(dac_addr_l,    daci       & 0x1 );                // select DAC
  digitalWrite(dac_addr_h,   (daci >> 1) & 0x1 );                // ----- || -----
//  if ( (dacpars[dac].bts & 0x02) == 0 ) llval = llval << 2;      // 18-bit
  spibuff[0] = (byte)((llval >> 16) & 0x0F) | 0x10;              // add 0x10 "write to DAC"
  spibuff[1] = (byte)((llval >> 8 ) & 0XFF);
  spibuff[2] = (byte)( llval        & 0XFF);
  digitalWrite(dac_SYNC, SYNC_LOW);                              // enable input register
  SPI.beginTransaction(SPISettings(spi_speed, spi_bit, spi_mode));  // ?
  SPI.transfer(spibuff,3);
  SPI.endTransaction();                                              // ?
  digitalWrite(dac_SYNC, SYNC_HIGH);                             // update input register
//  if ( (dacpars[dac].bts & 0x02) == 0 ) return llval >> 2;       // 18-bit
//     else                               return llval;            // 20-bit   
  return llval;             
}
*/

void serialEvent() {
  if (Serial.available()) {                         
    char inChar = (char)Serial.read();                            // get the new byte
    if (inputBuffer.length()==0) comtimer=millis();               // start timer for the first byte
    if ( (millis()-comtimer)>1000 ) inputBuffer="";               // reset input buffer after 1s delay 
    if (inChar == '\n') {                                         // end-of-command
      if ( inputBuffer == "IRQ" ) {
          IRQready = true;                                        // set interrupt flag
      } else {
          stringComplete = true;                                  // set complete command flug
          inputString=inputBuffer;                                // make a copy
      }    
      inputBuffer="";                                             // clear input buffer
    } else { inputBuffer += inChar; }                             // add char to the input Buffer:
    if (inputBuffer.length()==255) inputBuffer="";                // discard to prevent buffer overflow
  }
}

/*
 interrupt by the analog comparator
*/
void comp_irq() {
    cmp_flag = true;
}
