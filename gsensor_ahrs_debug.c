
#include <stm32f4xx_usart.h>

#include <stdlib.h>
#include <math.h>

#include <stdio.h>

#include "gsensor_ahrs.h"
#include "gsensor.h"
#include "sensor.h"

#include "csr8670_uart.h"

/*
  前進是ROLL-AXIS(Forward)(對y軸) 往旁邊偏是PITCH-AXIS(Size)(對z軸) 跳動是YAW-AXIS(Vertical)(對x軸)
  
  roll、pitch、yaw
  一個剛體有六個自由度：[x,y,y,θ,φ,ψ]
  其中的三個轉動自由度稱為：roll、pitch、yaw。常用的是後兩個，pitch又稱為tilt“傾斜”，
  yaw又被稱為pan“搖擺”。
  
  自由度—DOF
  A degree of freedom(DOF) refers to one of the set of independent position variabes, 
  with respect to a frame of reference, necessary to specify an object’s position within the world.
  自由度是指一組獨立的位置變數中的一個，它可以確定一個物體相對於一個參照系在空間中的位置。
  frame of reference—參照系
  機器人中的關節一般用“自由度”來表示。
  
  他們都是加速度 m(米)/s^2(秒平方)
  
  步數判斷就根據加速度變化 設一threshold來判斷之
*/


extern GSensor_LIS2DH12_t Axes_Data;

//===iNEMO_AHRS.c===================================//

#define iNEMO_CNT     (pSensorData->m_nCount)
#define iNEMO_DT      (pSensorData->m_fDeltaTime)

#define iNEMO_AccX    (pSensorData->m_fAcc[0])
#define iNEMO_AccY    (pSensorData->m_fAcc[1])
#define iNEMO_AccZ    (pSensorData->m_fAcc[2])
#define iNEMO_GyroX   (pSensorData->m_fGyro[0])
#define iNEMO_GyroY   (pSensorData->m_fGyro[1])
#define iNEMO_GyroZ   (pSensorData->m_fGyro[2])
#define iNEMO_MagX    (pSensorData->m_fMag[0])
#define iNEMO_MagY    (pSensorData->m_fMag[1])
#define iNEMO_MagZ    (pSensorData->m_fMag[2])
#define iNEMO_VarAcc  (pSensorData->m_fVarAcc)
#define iNEMO_VarMag  (pSensorData->m_fVarMag)
#define iNEMO_VarGyro (pSensorData->m_fVarGyro)
#define iNEMO_AccRef  (pSensorData->m_fAccRef) //X,Y,Z向量
#define iNEMO_MagRef  (pSensorData->m_fMagRef)

//==================================================//
#define LengthMax      4000
#define LengthStandard 3000
#define LengthMin      2000
//==================================================//

void float2String(char * DestnaStr ,double Data);

//==================================================//

//P = F * P F' + Q
iNEMO_fMATRIX_TYPE* pMat_P;    //Error Covariance Matrix  
iNEMO_fMATRIX_TYPE* pMat_Fsys; // The State Matrix (Jacobian)
iNEMO_fMATRIX_TYPE* pMat_Q;    //Covariance Matrix of the Noise Process
iNEMO_fMATRIX_TYPE* pMat_Pnew; //Error Covariance Matrix (updated)

iNEMO_fMATRIX_TYPE* pMat_Ka;   //Kalman Gain Matrix of Accelerometer
iNEMO_fMATRIX_TYPE* pMat_Km;   //Kalman Gain Matrix of Magnetometer

iNEMO_fMATRIX_TYPE* pMat_Hja;  //Jacobian Matrix of Accelerometer
iNEMO_fMATRIX_TYPE* pMat_Hjm;  //Jacobian Matrix of Magnetometer

iNEMO_fMATRIX_TYPE* pMat_Ra;   //The Accelerometer Measurement Noise
iNEMO_fMATRIX_TYPE* pMat_Rm;   //The Magnetometer Measurement Noise

float fSV[7];  //State Variable Vector  [4,5,6]=陀螺儀X,Y,Z
float fHalfDt; //The half of the DT

//X,Y,Z軸的平方相加開根號之值 用來校正用
float fAccRef; //Acceleration reference. The norm of this vector will be used for the accereation correction
float fMagRef;

//===AHRS_Continuous.c==============================//
float fAccXYZ[3]={0,0,0}, fMagXYZ[3]={0,0,0}; //Acceleration and Magnetic field values
float fGyroXYZ[3]={0,0,0}; //Gyroscope data

//AHRS filter input data structure
iNEMO_SENSORDATA xSensorData;
iNEMO_EULER_ANGLES xEulerAngles={0};
iNEMO_QUAT  xQuat={0};
uint8_t iter=0;

float t=0;
//==================================================//
/*
float AVGx, AVGz;
float xcont[7]={0, 0, 0, 0, 0, 0, 0}; 
float zcont[7]={0, 0, 0, 0, 0, 0, 0};
uint8_t iterj=0;
*/

float LengthAvg=0;

float step=0;

float Length=0; //矢量長度
uint8_t LengthState=0, lastLengthState=0; // state=0 :standard Lengh; state=1 :max/min Length
//==================================================//

//===iNEMO_AHRS_MemMan_1.c==========================//
void *iNEMO_Malloc(size_t size)
{ 
  return malloc(size);
}

void iNEMO_Free(void *p)
{
  free(p);
}
//==================================================//

//===iNEMO_AHRS.c===================================//
void iNEMO_AHRS_Init(iNEMO_SENSORDATA*    pSensorData,
                     iNEMO_EULER_ANGLES*  pAngle,
                     iNEMO_QUAT*          pQuat)
{

  /* Initialize all the EKF Matrix to Zero */
  
  pMat_P = iNEMO_fMatCreateZero(7,7);     //=iNEMO_fMatCreateInit(7, 7, NULL matrix), iNEMO_fMatCreateInit(Row, Col, Type)
  pMat_Fsys = iNEMO_fMatCreateUnit(7,7);  //=iNEMO_fMatCreateInit(7, 7, 全1的matrix)
  pMat_Q = iNEMO_fMatCreateZero(7,7);
  pMat_Pnew = iNEMO_fMatCreateZero(7,7);
  pMat_Ra = iNEMO_fMatCreateZero(3,3);
  pMat_Rm = iNEMO_fMatCreateZero(3,3);
  pMat_Ka = iNEMO_fMatCreateZero(7,3);
  pMat_Km = iNEMO_fMatCreateZero(7,3);
 
  /* Assign specific values within the matrixes */
  for(short int i=0; i < 7; ++i)
  {
    iNEMO_MatData(pMat_P)[i][i] = 1.0e-1f;  //對角線給值 =(1e-1)f=1/10=0.1
    iNEMO_MatData(pMat_Q)[i][i] = 5.0e-7f;  //=(5e-7)f = 5/10000000=0.0000005
  }


  for(short int i=0; i<3; i++)
  {
    iNEMO_MatData(pMat_Ra)[i][i] = iNEMO_VarAcc;  //對角線給值: G-Sensor的誤差值
    
    iNEMO_MatData(pMat_Rm)[i][i] = iNEMO_VarMag;
  }
  
  pMat_Hja = iNEMO_fMatCreateZero(3,7);
  pMat_Hjm = iNEMO_fMatCreateZero(3,7);


  /* State Vector initialization... */
  /* ...Option 1: {1, 0, 0, 0, 0, 0, 0} */
  fSV[0] = 1;
  for(short int i=1; i < 7; ++i)
    fSV[i] = 0;
  
  /* Sensor Data Initialization */
  iNEMO_CNT = 0;
  
  iNEMO_AccX = 0.0f;
  iNEMO_AccY = 0.0f;
  iNEMO_AccZ = 0.0f;
  iNEMO_GyroX = 0.0f;
  iNEMO_GyroY = 0.0f;
  iNEMO_GyroZ = 0.0f;
  iNEMO_MagX = 0.0f;
  iNEMO_MagY = 0.0f;
  iNEMO_MagZ = 0.0f;
  
  /* Angle initialization... */
  /* ...Option 1 (equivalent to O1 of State Vector) */
  iNEMO_Roll(pAngle) = 0.0f;
  iNEMO_Pitch(pAngle) = 0.0f;
  iNEMO_Yaw(pAngle) = 0.0f;
  
  fHalfDt=iNEMO_DT/2;
  
  //(AccRef[0]^2+AccRef[1]^2+AccRef[2]^2)開根號 = 新的fAccRef  (AccRef[0..2] = X,Y,Z軸)
  fAccRef=0.0;
  fMagRef=0.0;
  for(short int i=0; i<3; i++)
  {
    fAccRef+=(iNEMO_AccRef[i]*iNEMO_AccRef[i]); 
    fMagRef+=(iNEMO_MagRef[i]*iNEMO_MagRef[i]);
  }
  
  fAccRef=sqrt(fAccRef);
  fMagRef=sqrt(fMagRef);
  
}

void iNEMO_AHRS_DeInit(iNEMO_SENSORDATA*    pSensorData,
                       iNEMO_EULER_ANGLES*  pAngle,
                       iNEMO_QUAT*          pQuat)
{
  /* Delete the basic EKF Matrix */
  iNEMO_fMatFree(pMat_P);
  iNEMO_fMatFree(pMat_Fsys);
  iNEMO_fMatFree(pMat_Q);
  
  iNEMO_fMatFree(pMat_Ra);
  iNEMO_fMatFree(pMat_Rm);
  iNEMO_fMatFree(pMat_Ka);
  iNEMO_fMatFree(pMat_Km);
  
  iNEMO_fMatFree(pMat_Hja);
  iNEMO_fMatFree(pMat_Hjm);
}

void iNEMO_AHRS_Update(iNEMO_SENSORDATA*    pSensorData,
                       iNEMO_EULER_ANGLES*  pAngle,
                       iNEMO_QUAT*          pQuat)
{
  int i;
  float fGyroX, fGyroY, fGyroZ;
  iNEMO_QUAT fSVnorm; //fSV的四元數: typedef float iNEMO_QUAT[4];
  float fTempNorm;
  float fH[4];

  /* Self increment of the Counter of the Update Sessions */
  iNEMO_CNT++;
  
  /* calculate Gyro Values */
  fGyroX = (iNEMO_GyroX - fSV[4]) * fHalfDt;
  fGyroY = (iNEMO_GyroY - fSV[5]) * fHalfDt;
  fGyroZ = (iNEMO_GyroZ - fSV[6]) * fHalfDt;
  
  /*state transition matrix */
  iNEMO_MatData(pMat_Fsys)[0][1] = -fGyroX;
  iNEMO_MatData(pMat_Fsys)[0][2] = -fGyroY;
  iNEMO_MatData(pMat_Fsys)[0][3] = -fGyroZ;
  iNEMO_MatData(pMat_Fsys)[1][0] =  fGyroX;
  iNEMO_MatData(pMat_Fsys)[1][2] =  fGyroZ;
  iNEMO_MatData(pMat_Fsys)[1][3] = -fGyroY;
  iNEMO_MatData(pMat_Fsys)[2][0] =  fGyroY;
  iNEMO_MatData(pMat_Fsys)[2][1] = -fGyroZ;
  iNEMO_MatData(pMat_Fsys)[2][3] =  fGyroX;
  iNEMO_MatData(pMat_Fsys)[3][0] =  fGyroZ;
  iNEMO_MatData(pMat_Fsys)[3][1] =  fGyroY;
  iNEMO_MatData(pMat_Fsys)[3][2] = -fGyroX;
  
  iNEMO_MatData(pMat_Fsys)[0][4] =  fSV[1] * fHalfDt;
  iNEMO_MatData(pMat_Fsys)[0][5] =  fSV[2] * fHalfDt;
  iNEMO_MatData(pMat_Fsys)[0][6] =  fSV[3] * fHalfDt;
  iNEMO_MatData(pMat_Fsys)[1][4] = -fSV[0] * fHalfDt;
  iNEMO_MatData(pMat_Fsys)[1][5] =  fSV[3] * fHalfDt;
  iNEMO_MatData(pMat_Fsys)[1][6] = -fSV[2] * fHalfDt;
  iNEMO_MatData(pMat_Fsys)[2][4] = -fSV[3] * fHalfDt;
  iNEMO_MatData(pMat_Fsys)[2][5] = -fSV[0] * fHalfDt;
  iNEMO_MatData(pMat_Fsys)[2][6] =  fSV[1] * fHalfDt;
  iNEMO_MatData(pMat_Fsys)[3][4] =  fSV[2] * fHalfDt;
  iNEMO_MatData(pMat_Fsys)[3][5] = -fSV[1] * fHalfDt;
  iNEMO_MatData(pMat_Fsys)[3][6] = -fSV[0] * fHalfDt;
  
  /*
  ***********************************************************************
  Extended Kalman Filter: Prediction Step
  ***********************************************************************
  */
  
  /* Update Quaternion with the new gyroscope measurements */
  fSVnorm[0] = fSV[0] - (fGyroX * fSV[1]) - (fGyroY * fSV[2]) - (fGyroZ * fSV[3]);
  fSVnorm[1] = fSV[1] + (fGyroX * fSV[0]) - (fGyroY * fSV[3]) + (fGyroZ * fSV[2]);
  fSVnorm[2] = fSV[2] + (fGyroX * fSV[3]) + (fGyroY * fSV[0]) - (fGyroZ * fSV[1]);
  fSVnorm[3] = fSV[3] - (fGyroX * fSV[2]) + (fGyroY * fSV[1]) + (fGyroZ * fSV[0]);
  
  for (i=0; i < 4; ++i)
    fSV[i] = fSVnorm[i];
  
  /*  P = F * P F' + Q  */
  pMat_Pnew = iNEMO_PropagateP(pMat_P, pMat_Fsys, pMat_Q, pMat_Pnew);
  
  /* Copy the new P in P */
  iNEMO_fMatCopy(pMat_Pnew, pMat_P);
  
  /*
  ***********************************************************************
  Extended Kalman Filter: Correction Step for tilting(v.傾斜,擺動) m_angle
  ***********************************************************************
  Perform the Correction any two updates 每兩次update就做一次校正
   */
  
  if(iNEMO_CNT % 2 == 0)
  {
    /* Assing the non-linear function h, it relates the State Variables with
    the measurements */
    fH[0] = -2.0*fAccRef * (fSV[1] * fSV[3] - fSV[0] * fSV[2]);
    fH[1] = -2.0*fAccRef * (fSV[0] * fSV[1] + fSV[2] * fSV[3]);
    fH[2] = -fAccRef  * (fSV[0] * fSV[0] - fSV[1] * fSV[1] - fSV[2] * fSV[2] + fSV[3] * fSV[3]);
  
    /* Its Jacobian matrix */
    iNEMO_MatData(pMat_Hja)[0][0] =  2.0*fAccRef * fSV[2];
    iNEMO_MatData(pMat_Hja)[0][1] = -2.0*fAccRef * fSV[3];
    iNEMO_MatData(pMat_Hja)[0][2] =  2.0*fAccRef * fSV[0];
    iNEMO_MatData(pMat_Hja)[0][3] = -2.0*fAccRef * fSV[1];
  
    iNEMO_MatData(pMat_Hja)[1][0] = -2.0*fAccRef * fSV[1];
    iNEMO_MatData(pMat_Hja)[1][1] = -2.0*fAccRef * fSV[0];
    iNEMO_MatData(pMat_Hja)[1][2] = -2.0*fAccRef * fSV[3];
    iNEMO_MatData(pMat_Hja)[1][3] = -2.0*fAccRef * fSV[2];
  
    iNEMO_MatData(pMat_Hja)[2][0] = -2.0*fAccRef * fSV[0];
    iNEMO_MatData(pMat_Hja)[2][1] =  2.0*fAccRef * fSV[1];
    iNEMO_MatData(pMat_Hja)[2][2] =  2.0*fAccRef * fSV[2];
    iNEMO_MatData(pMat_Hja)[2][3] = -2.0*fAccRef * fSV[3];
  
    /* Ka = P * Hja' * (Hja * P * Hja' + Ra)^-1  */
    pMat_Ka = iNEMO_CalculateK(pMat_P, pMat_Hja, pMat_Ra, pMat_Ka);
  
    /* Update State Vector */
    for (i=0; i < 7; ++i)
      fSV[i] += iNEMO_MatData(pMat_Ka)[i][0] * (iNEMO_AccX - fH[0]) +
        iNEMO_MatData(pMat_Ka)[i][1] * (iNEMO_AccY - fH[1]) +
          iNEMO_MatData(pMat_Ka)[i][2] * (iNEMO_AccZ - fH[2]);
  
    /* P = (I - Ka*Hja)*P */
    pMat_Pnew = iNEMO_UpdateP(pMat_P, pMat_Ka, pMat_Hja, pMat_Pnew);
    /* Copy the new P in P */
    iNEMO_fMatCopy(pMat_Pnew, pMat_P);
  } /* endif CNT /2 */

  /*
  ***********************************************************************
  Extended Kalman Filter: Correction Step for heading(航向) m_angle
  ***********************************************************************
  Perform the Correction any two updates
  */
  if(iNEMO_CNT % 2  == 0)
  {
    /* Yaw Correction */
  
    /* Assing the non-linear function h, it relates the State Variables with
    the measurements */
    fH[0] =         fMagRef * (fSV[0]*fSV[0] + fSV[1]*fSV[1] - fSV[2]*fSV[2] - fSV[3]*fSV[3]);
    fH[1] =  2.0f * fMagRef * (fSV[1]*fSV[2] - fSV[3]*fSV[0]);
    fH[2] =  2.0f * fMagRef * (fSV[1]*fSV[3] + fSV[2]*fSV[0]);
  
    /* Calculate the relative Jacobian Matrix */
    iNEMO_MatData(pMat_Hjm)[0][0] =   2.0f * fMagRef * fSV[0];
    iNEMO_MatData(pMat_Hjm)[0][1] =   2.0f * fMagRef * fSV[1];
    iNEMO_MatData(pMat_Hjm)[0][2] = - 2.0f * fMagRef * fSV[2];
    iNEMO_MatData(pMat_Hjm)[0][3] = - 2.0f * fMagRef * fSV[3];
    iNEMO_MatData(pMat_Hjm)[1][0] = - 2.0f * fMagRef * fSV[3];
    iNEMO_MatData(pMat_Hjm)[1][1] =   2.0f * fMagRef * fSV[2];
    iNEMO_MatData(pMat_Hjm)[1][2] =   2.0f * fMagRef * fSV[1];
    iNEMO_MatData(pMat_Hjm)[1][3] = - 2.0f * fMagRef * fSV[0];
    iNEMO_MatData(pMat_Hjm)[2][0] =   2.0f * fMagRef * fSV[2];
    iNEMO_MatData(pMat_Hjm)[2][1] =   2.0f * fMagRef * fSV[3];
    iNEMO_MatData(pMat_Hjm)[2][2] =   2.0f * fMagRef * fSV[0];
    iNEMO_MatData(pMat_Hjm)[2][3] =   2.0f * fMagRef * fSV[1];
  
    /* Km = P * Hjm' * (Hjm * P * Hjm' + Rm)^-1  */
    pMat_Km = iNEMO_CalculateK(pMat_P, pMat_Hjm, pMat_Rm, pMat_Km);
  
    /* State Vector Update */
    for (i=0; i < 7; ++i)
      fSV[i] += iNEMO_MatData(pMat_Km)[i][0] * (iNEMO_MagX - fH[0]) +
        iNEMO_MatData(pMat_Km)[i][1] * (iNEMO_MagY - fH[1]) +
          iNEMO_MatData(pMat_Km)[i][2] * (iNEMO_MagZ - fH[2]);
  
    /* Update of Error Covariance Matrix */
    /* P = (I - Km*Hjm*P */
    pMat_Pnew = iNEMO_UpdateP(pMat_P, pMat_Km, pMat_Hjm, pMat_Pnew);
    /* Copy the new P in P */
    iNEMO_fMatCopy(pMat_Pnew, pMat_P);
  
  } /* End of Yaw Correction */
  
  /* Rescale quaternion to have module 1 */
  fTempNorm = 1.0f / sqrt(fSV[0]*fSV[0] + fSV[1]*fSV[1] + fSV[2]*fSV[2] + fSV[3]*fSV[3]);
  for(i = 0; i < 4; ++i)
  {
    fSV[i] = fTempNorm * fSV[i];
    pQuat[0][i]=fSV[i];
  }
  
  /* Update the RPY angles */
  iNEMO_Roll(pAngle)  = atan2(2.0f*(fSV[0]*fSV[1] + fSV[2]*fSV[3]), 1.0f - 2.0f*(fSV[1]*fSV[1] + fSV[2]*fSV[2]));
  
  iNEMO_Pitch(pAngle) = asin( -2.0f*(fSV[1]*fSV[3] - fSV[0]*fSV[2]));
  
  iNEMO_Yaw(pAngle)   = atan2(2.0f*(fSV[1]*fSV[2] + fSV[0]*fSV[3]), (1.0f - 2.0f*(fSV[2]*fSV[2] + fSV[3]*fSV[3])));
  
  
}

//==================================================//
//===iNEMO_EKF.c====================================//
iNEMO_fMATRIX_TYPE* iNEMO_PropagateP(iNEMO_fMATRIX_TYPE* pPoldMat, 
                                     iNEMO_fMATRIX_TYPE* pStateMat, 
                                     iNEMO_fMATRIX_TYPE* pQMat, 
                                     iNEMO_fMATRIX_TYPE* pPnewMat)
{
  /* Internal Variables */
  iNEMO_fMATRIX_TYPE* pTmp;

  pTmp = iNEMO_fMatCreate(iNEMO_MatRow(pPoldMat), iNEMO_MatCol(pPoldMat));
  
  /* First Step: pTmp = pStateMat * pPoldMat */
  pTmp = iNEMO_fMatMulMat(pStateMat, pPoldMat, pTmp);

  /* Second Step: pPnewMat = pTmp * pStateMat' */
  pPnewMat = iNEMO_fMatMulMatMT(pTmp, pStateMat, pPnewMat);

  /* Third Step: pPnewMat += Q */
  pPnewMat = iNEMO_fMatAdd(pPnewMat, pQMat, pPnewMat);
  
  /* Delete the Temporary Matrix before to exit */
  iNEMO_fMatFree(pTmp);
  
  return pPnewMat;  
}

iNEMO_fMATRIX_TYPE* iNEMO_CalculateK(iNEMO_fMATRIX_TYPE* pPMat,
                                     iNEMO_fMATRIX_TYPE* pHJMat,
                                     iNEMO_fMATRIX_TYPE* pRMat,
                                     iNEMO_fMATRIX_TYPE* pKMat)
{
  /* Internal Variables */
  iNEMO_fMATRIX_TYPE* pTmp1;
  iNEMO_fMATRIX_TYPE* pTmp2;
  iNEMO_fMATRIX_TYPE* pRInvMat;
  
  pTmp1 = iNEMO_fMatCreate(iNEMO_MatRow(pPMat), iNEMO_MatRow(pHJMat));
  pTmp2 = iNEMO_fMatCreate(iNEMO_MatRow(pHJMat), iNEMO_MatCol(pTmp1));  
  
  /* First Step: pTmp1 = pPMat * pHJMat' */
  pTmp1 = iNEMO_fMatMulMatMT(pPMat, pHJMat, pTmp1);
  
  /* Second Step: pTmp2 = pHJMat * (pPMat * pHJMat') = pHJMat * pTmp1 */
  pTmp2 = iNEMO_fMatMulMat(pHJMat, pTmp1, pTmp2);
  
  /* Third Step: pTmp2 = pTmp2 + pRMat */
  pTmp2 = iNEMO_fMatAdd(pTmp2, pRMat, pTmp2);
  
  /* Fourth Step: pRInvMat = (pTmp2)^-1 */
  pRInvMat = iNEMO_fMatCreate(iNEMO_MatRow(pTmp2), iNEMO_MatCol(pTmp2));
  
  pRInvMat = iNEMO_fMatInv(pTmp2, pRInvMat);
  
  /* Fifth Step: pKMat = pTmp1 * pRInvMat */
  pKMat = iNEMO_fMatMulMat(pTmp1, pRInvMat, pKMat);
  
  /* Delete the Temporary Matrix before to exit */
  iNEMO_fMatFree(pTmp1);
  iNEMO_fMatFree(pTmp2);  
  iNEMO_fMatFree(pRInvMat);  
  
  return pKMat;
  
}

iNEMO_fMATRIX_TYPE* iNEMO_UpdateP(iNEMO_fMATRIX_TYPE* pPoldMat,
                                  iNEMO_fMATRIX_TYPE* pKMat,
                                  iNEMO_fMATRIX_TYPE* pHJMat, 
                                  iNEMO_fMATRIX_TYPE* pPnewMat)
{
  iNEMO_fMATRIX_TYPE* pTmp;
  iNEMO_fMATRIX_TYPE* pI;
   
  pTmp = iNEMO_fMatCreateZero(iNEMO_MatRow(pPoldMat), iNEMO_MatCol(pPoldMat));
  pI = iNEMO_fMatCreateUnit(iNEMO_MatRow(pPoldMat), iNEMO_MatCol(pPoldMat));

  // First Step: pTmp = (pKMat * pHJMat)
  pTmp = iNEMO_fMatMulMat(pKMat, pHJMat, pTmp); 

  // Second Step: pTmp = I - pTmp
  pTmp = iNEMO_fMatSub(pI, pTmp, pTmp);
  
  // Third Step: pPnewMat = pTmp * pPoldMat
  pPnewMat = iNEMO_fMatMulMat(pTmp, pPoldMat, pPnewMat);
  
  iNEMO_fMatFree(pTmp); 
  iNEMO_fMatFree(pI);
  
  return (pPnewMat);  
}
//===iNEMO_math.c===================================//
iNEMO_fMATRIX_TYPE *iNEMO_fMatCreate(int nRow, int nCol)
{
  int i;

  iNEMO_fMATRIX_TYPE *pTmp = (iNEMO_fMATRIX_TYPE*)iNEMO_Malloc(sizeof(iNEMO_fMATRIX_TYPE));

  if(pTmp != NULL)
  {
    // Check on null values
    if(nRow==0 || nCol==0)
      return(NULL);

    pTmp->m_nRow = nRow;
    pTmp->m_nCol = nCol;

    // Allocate memory for data
    iNEMO_MatData(pTmp) = iNEMO_Malloc(nRow * sizeof(float*));
    for(i=0; i < nRow; ++i)
      iNEMO_MatData(pTmp)[i] = iNEMO_Malloc(nCol * sizeof(float));
  }

  return pTmp;
}

iNEMO_sMATRIX_TYPE *iNEMO_sMatCreate(int nRow, int nCol)
{
  int i;
  iNEMO_sMATRIX_TYPE *pTmp = (iNEMO_sMATRIX_TYPE*)iNEMO_Malloc(sizeof(iNEMO_sMATRIX_TYPE));
  
  if(pTmp != NULL)
  {
    // Check on null values
    if(nRow==0 || nCol==0)
      return(NULL);
  
    pTmp->m_nRow = nRow;
    pTmp->m_nCol = nCol;
  
    // Allocate memory for data
    iNEMO_MatData(pTmp) = iNEMO_Malloc(nRow * sizeof(short int*));
    for(i=0; i < nRow; ++i)
      iNEMO_MatData(pTmp)[i] = iNEMO_Malloc(nCol * sizeof(short int));
  }
  
  return pTmp;
}

iNEMO_fMATRIX_TYPE *iNEMO_fMatCreateInit(int nRow, int nCol, int nType)
{
  iNEMO_fMATRIX_TYPE* pTmp;
  int i, j;
  
  if(nRow==0 || nCol==0)
    return(NULL);
  
  if((pTmp = iNEMO_fMatCreate( nRow, nCol )) != NULL)
  {
    switch(nType)
    {
      case iNEMO_ZERO_MATRIX:
      case iNEMO_IDEN_MATRIX:
      case iNEMO_ONES_MATRIX:
        for(i=0; i < nRow; i++)
        {
          for(j=0; j < nCol; j++)
          {
            /* Fills of ones */
            if(nType == iNEMO_ONES_MATRIX)
            {
              iNEMO_MatData(pTmp)[i][j] = 1.0f;
              continue;
            }
            /* Fills as Identity Matrix */
            if(nType == iNEMO_IDEN_MATRIX)
            {
              if(i==j)
              {
                iNEMO_MatData(pTmp)[i][j] = 1.0f;
                continue;
              }
            }
            /* Fills of zeros */
            iNEMO_MatData(pTmp)[i][j] = 0.0f;
          }
        }
        break;
    }
    return (pTmp);
  }
  else
    return (NULL);
}

iNEMO_fMATRIX_TYPE *iNEMO_fMatFill(iNEMO_fMATRIX_TYPE *pMat, float fValue)
{
  int i, j;

  for(i=0; i<iNEMO_MatRow(pMat); i++)
    for(j=0; j<iNEMO_MatCol(pMat); j++)
      iNEMO_MatData(pMat)[i][j] = (float) (fValue);
  return (pMat);
}

int iNEMO_fMatFree(iNEMO_fMATRIX_TYPE *pMat)
{
  int i;

  if(pMat == NULL)
    return (0);

  for(i=0; i<iNEMO_MatRow(pMat); i++)
  {
    iNEMO_Free( iNEMO_MatData(pMat)[i] );
  }
  iNEMO_Free(iNEMO_MatData(pMat));
  
  iNEMO_Free( pMat );
  
  return (1);
}

int iNEMO_sMatFree(iNEMO_sMATRIX_TYPE *pMat)
{
  int i;
  
  if(pMat == NULL)
    return (0);

  for(i=0; i<iNEMO_MatRow(pMat); i++)
  {
    iNEMO_Free( iNEMO_MatData(pMat)[i] );
  }
  iNEMO_Free(iNEMO_MatData(pMat));
  iNEMO_Free( pMat );
  return (1);
}

iNEMO_fMATRIX_TYPE* iNEMO_fMatCopy(iNEMO_fMATRIX_TYPE* pSource, iNEMO_fMATRIX_TYPE* pDest)
{
  int  i, j;
  
  // Check the dimensions
  if( iNEMO_MatRow(pSource) != iNEMO_MatRow(pDest) || iNEMO_MatCol(pSource) != iNEMO_MatCol(pDest) )
    return NULL;
  else
  {
    for (i=0; i < iNEMO_MatRow(pSource); i++)
      for (j=0; j < iNEMO_MatCol(pSource); j++)
        iNEMO_MatData(pDest)[i][j] = iNEMO_MatData(pSource)[i][j];
    return(pDest);
  }
}

iNEMO_fMATRIX_TYPE* iNEMO_fMatAdd(iNEMO_fMATRIX_TYPE* pTerm1,
                                  iNEMO_fMATRIX_TYPE* pTerm2,
                                  iNEMO_fMATRIX_TYPE* pAdd)
{
  int	i, j;
  
  /* Check if dimensions are wrong */
  if( iNEMO_MatRow(pAdd) != iNEMO_MatRow(pTerm1) || iNEMO_MatCol(pAdd) != iNEMO_MatCol(pTerm2) )
    return NULL;
  else
  {
	  for(i=0; i < iNEMO_MatRow(pTerm1); ++i)
      for(j=0; j < iNEMO_MatCol(pTerm1); ++j)
        iNEMO_MatData(pAdd)[i][j] = iNEMO_MatData(pTerm1)[i][j] + iNEMO_MatData(pTerm2)[i][j];
  }
  return(pAdd);
}

iNEMO_fMATRIX_TYPE* iNEMO_fMatSub(iNEMO_fMATRIX_TYPE* pTerm1,
                                  iNEMO_fMATRIX_TYPE* pTerm2,
                                  iNEMO_fMATRIX_TYPE* pSub)
{
  int	i, j;

  /* Check if dimensions are wrong */
  if( iNEMO_MatRow(pTerm1) != iNEMO_MatRow(pSub) || iNEMO_MatCol(pTerm1) != iNEMO_MatCol(pSub) )
    return NULL;
  else
  {
    for (i=0; i < iNEMO_MatRow(pTerm1); ++i)
      for (j=0; j < iNEMO_MatCol(pTerm1); ++j)
        iNEMO_MatData(pSub)[i][j] = iNEMO_MatData(pTerm1)[i][j] - iNEMO_MatData(pTerm2)[i][j];
  }
  return(pSub);
}

// 橫的 乘以 直的
iNEMO_fMATRIX_TYPE* iNEMO_fMatMulMat(iNEMO_fMATRIX_TYPE* pTerm1,
                                     iNEMO_fMATRIX_TYPE* pTerm2,
                                     iNEMO_fMATRIX_TYPE* pMul)
{
  int	i, j, k;

  // if dimensions are wrong
  if( iNEMO_MatRow(pMul) != iNEMO_MatRow(pTerm1) || iNEMO_MatCol(pMul) != iNEMO_MatCol(pTerm2) )
    return NULL;
  else
  {
    for(i=0; i < iNEMO_MatRow(pTerm1); i++)
      for(j=0; j < iNEMO_MatCol(pTerm2); j++)
        for(k=0, iNEMO_MatData(pMul)[i][j]=0.0; k < iNEMO_MatCol(pTerm1); k++)
          iNEMO_MatData(pMul)[i][j] += iNEMO_MatData(pTerm1)[i][k] * iNEMO_MatData(pTerm2)[k][j];
  }

  return(pMul);
}

// 橫的 乘以 橫的
iNEMO_fMATRIX_TYPE* iNEMO_fMatMulMatMT(iNEMO_fMATRIX_TYPE* pTerm1,
                                       iNEMO_fMATRIX_TYPE* pTerm2,
                                       iNEMO_fMATRIX_TYPE* pMul)
{
  int i, j, k;

  // if dimensions are wrong
  if( iNEMO_MatRow(pMul) != iNEMO_MatRow(pTerm1) || iNEMO_MatCol(pMul) != iNEMO_MatRow(pTerm2) )
    return NULL;
  else
  {
    for (i=0; i < iNEMO_MatRow(pTerm1); i++)
      for (j=0; j < iNEMO_MatRow(pTerm2); j++)
        for (k=0, iNEMO_MatData(pMul)[i][j]=0.0; k < iNEMO_MatCol(pTerm1); k++)
          iNEMO_MatData(pMul)[i][j] += iNEMO_MatData(pTerm1)[i][k] * iNEMO_MatData(pTerm2)[j][k];
  }
  return(pMul);
}

iNEMO_fMATRIX_TYPE* iNEMO_fMatInv(iNEMO_fMATRIX_TYPE* pSource,
                                  iNEMO_fMATRIX_TYPE* pDest)
{
  iNEMO_fMATRIX_TYPE* A;
  iNEMO_fMATRIX_TYPE* B;
  iNEMO_sMATRIX_TYPE* P;
  int i, nCol, nRow;
  
  nCol = iNEMO_MatCol(pSource);
  nRow = iNEMO_MatRow(pSource);
  
  if(nCol != nRow)
    /* The matrix is not Square */
    return (NULL);
  
  A = iNEMO_fMatCreate(nRow, nCol);
  if(A == NULL)
    return (NULL);
  
  B = iNEMO_fMatCreate(nRow, nCol);
  if(B == NULL)
    return (NULL);
  
  /* P is a vector matrix */
  P = iNEMO_sMatCreate(nRow, 1);
  if(P == NULL)
    return (NULL);
  
  /* It is to avoid to modify pSource Matrix */
  iNEMO_fMatCopy(pSource, A);
  
  /* LU Decomposition and check for Singular Matrix */
  if (iNEMO_MatLUP(A, P) == -1)
  {
    iNEMO_fMatFree(A);
    iNEMO_fMatFree(B);
    iNEMO_sMatFree(P);
  
    return (NULL);
  }
  
  for (i=0; i<nCol; ++i)
  {
    iNEMO_fMatFill(B, 0.0f);
  
    iNEMO_MatData(B)[i][0] = 1.0f;
    iNEMO_MatBackSubs(A, B, P, pDest, i);
  }
  iNEMO_fMatFree(A);
  iNEMO_fMatFree(B);
  iNEMO_sMatFree(P);
  
  if(pDest == NULL)
  {
    return(NULL);
  }
  else
  {
    return (pDest);
  }
}

int iNEMO_MatLUP(iNEMO_fMATRIX_TYPE* pSourceDestLU, iNEMO_sMATRIX_TYPE* pPerm)
{
  int	i, j, k, iC;
  int	iMax;
  int	retNumPerm = 0;
  short int sTmp;
  float fP1, fP2; /* Pivot Variables */
  
  iC = iNEMO_MatCol(pSourceDestLU);
  
  for(i=0; i < iC; ++i)
    iNEMO_MatData(pPerm)[i][0] = (short int) (i);
  
  /* Partial Pivoting */
  for (k=0; k < iC; ++k)
  {
    for (i=k, iMax=k, fP1=0.0f; i < iC; ++i)
    {
      /* Local ABS */
      if (iNEMO_MatData(pSourceDestLU)[iNEMO_MatData(pPerm)[i][0]][k] > 0)
        fP2 = iNEMO_MatData(pSourceDestLU)[iNEMO_MatData(pPerm)[i][0]][k];
      else
        fP2 = - iNEMO_MatData(pSourceDestLU)[iNEMO_MatData(pPerm)[i][0]][k];
      if(fP2 > fP1)
      {
        fP1 = fP2;
        iMax = i;
      }
    }
    /* Row exchange, update permutation vector */
    if(k != iMax)
    {
      retNumPerm++;
      sTmp = iNEMO_MatData(pPerm)[k][0];
      iNEMO_MatData(pPerm)[k][0] = iNEMO_MatData(pPerm)[iMax][0];
      iNEMO_MatData(pPerm)[iMax][0] = sTmp;
    }
    
    /* Suspected Singular Matrix */
    if(iNEMO_MatData(pSourceDestLU)[iNEMO_MatData(pPerm)[k][0]][k] == 0.0f)
      return (-1);
    
    for(i=k+1; i < iC; ++i)
    {
      /* Calculate Mat [i][j] */
      iNEMO_MatData(pSourceDestLU)[iNEMO_MatData(pPerm)[i][0]][k] =
        iNEMO_MatData(pSourceDestLU)[iNEMO_MatData(pPerm)[i][0]][k] /
        iNEMO_MatData(pSourceDestLU)[iNEMO_MatData(pPerm)[k][0]][k];
      
      /* Elimination */
      for(j=k+1; j < iC; ++j)
      {
        iNEMO_MatData(pSourceDestLU)[iNEMO_MatData(pPerm)[i][0]][j] -=
          iNEMO_MatData(pSourceDestLU)[iNEMO_MatData(pPerm)[i][0]][k] *
          iNEMO_MatData(pSourceDestLU)[iNEMO_MatData(pPerm)[k][0]][j];
      }
      
      
    }
  }
  return (retNumPerm);
}

iNEMO_fMATRIX_TYPE* iNEMO_MatBackSubs(iNEMO_fMATRIX_TYPE* pSourceLU,
                                      iNEMO_fMATRIX_TYPE* pSourceDestColumn,
                                      iNEMO_sMATRIX_TYPE* pPerm,
                                      iNEMO_fMATRIX_TYPE* pDest,
                                      int iResultCol)
{
  int i, j, k, iC;
  float fSum, fTmp;
  
  iC = iNEMO_MatCol(pSourceLU);
  
  for(k=0; k < iC; ++k)
    for(i=k+1; i < iC; ++i)
      iNEMO_MatData(pSourceDestColumn)[iNEMO_MatData(pPerm)[i][0]][0] -=
        iNEMO_MatData(pSourceLU)[iNEMO_MatData(pPerm)[i][0]][k] *
          iNEMO_MatData(pSourceDestColumn)[iNEMO_MatData(pPerm)[k][0]][0];
  
  iNEMO_MatData(pDest)[iC-1][iResultCol] = iNEMO_MatData(pSourceDestColumn)[iNEMO_MatData(pPerm)[iC-1][0]][0] /
    iNEMO_MatData(pSourceLU)[iNEMO_MatData(pPerm)[iC-1][0]][iC-1];
  
  
  for(k=iC-2; k >= 0; k--)
  {
    fSum = 0.0f;
    
    for (j=k+1; j < iC; ++j)
      fSum += iNEMO_MatData(pSourceLU)[iNEMO_MatData(pPerm)[k][0]][j] * iNEMO_MatData(pDest)[j][iResultCol];
    
    fTmp = iNEMO_MatData(pSourceDestColumn)[iNEMO_MatData(pPerm)[k][0]][0] - fSum;
    iNEMO_MatData(pDest)[k][iResultCol] = fTmp /
     iNEMO_MatData(pSourceLU)[iNEMO_MatData(pPerm)[k][0]][k];
  }
  
  return pDest;
  
}

float iNEMO_WrapAround(float fInput)
{
  if(fInput > PI)
    fInput -= PI2;
    
  if(fInput < -PI)
    fInput += PI2;
    
  return (fInput);
}
//==================================================//

//===AHRS_Continuous.c==============================//
void GSensorPedometer_Init(void)
{
  /* Filter references for Acceleration and Magnetic field */
  xSensorData.m_fAccRef[0]=0;
  xSensorData.m_fAccRef[1]=0;
  xSensorData.m_fAccRef[2]=-9.81f; //這邊可能要做G值轉換: 1G=9.81
      
  xSensorData.m_fMagRef[0]=0.37f;
  xSensorData.m_fMagRef[1]=0;
  xSensorData.m_fMagRef[2]=-0.25f;
  
  xSensorData.m_fDeltaTime=0.02f;
  
  xSensorData.m_fVarAcc=5.346e-6;
  xSensorData.m_fVarMag=5.346e-6;

  iNEMO_AHRS_Init(&xSensorData, &xEulerAngles, &xQuat);  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
void float2String(char * DestnaStr, double Data)
{
  signed int DataZ; //整數部分
  
  DataZ  = (int)Data; //提整數部分 

  if(!DataZ )  //考慮特殊情況
    *DestnaStr++ = '0';
  else
  {
    //負值考慮
    if(DataZ<0)
    {
      *DestnaStr = 0x2D;
      DestnaStr++;

      //變正值
      DataZ = 0 - DataZ;
      
      while(DataZ )
      {
        *DestnaStr = DataZ  % 10 + 48; //一次讀整數每一位
        DataZ  = DataZ  / 10;
        DestnaStr++;  //儲存下一個位置
      }
      
    }
    else
    {
      while(DataZ )
      {
        *DestnaStr = DataZ  % 10 + 48; //一次讀整數每一位
        DataZ  = DataZ  / 10;
        DestnaStr++;  //儲存下一個位置
      }
    }  
  }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int GSensorPedometer(void)
{
  char DestnaStrLen[10];
  char DestnaStrX[10];
  char DestnaStrY[10];
  char DestnaStrZ[10];
  char DestnaStrT[10];
  char DestnaStrStep[10];
  char DestnaStrLengthAvg[10];
  int i=0;

  //float temp[2];
  
  //float roll, pitch, yaw;
  /*
  for(i=0; i<10; i++)
  {
  	DestnaStrLen[i]=0;
  	DestnaStrX[i]=0;
  	DestnaStrY[i]=0;
  	DestnaStrZ[i]=0;
    DestnaStrT[i]=0;
    DestnaStrStep[i]=0; 
  }
  // Read data
  fAccXYZ[0] = toG((float)Axes_Data.X) * 1000; //mg單位
  fAccXYZ[1] = toG((float)Axes_Data.Y) * 1000;
  fAccXYZ[2] = toG((float)Axes_Data.Z) * 1000; 
  */
  
  //fAccXYZ[0] = (float)(Axes_Data.X); //2 g full scale [LSB/mg]=>LSM_Acc_Sensitivity_2g=1.0  (Axes_Data->X)單位LSB / LSM_Acc_Sensitivity_2g = X 單位mg
  //fAccXYZ[1] = (float)(Axes_Data.Y);
  //fAccXYZ[2] = (float)(Axes_Data.Z);
  /*
  temp[0]= -fGyroXYZ[1];
  temp[1] = fGyroXYZ[0];
  fGyroXYZ[0] = temp[0];
  fGyroXYZ[1] = temp[1];  
  
  {
    xSensorData.m_fAcc[0]=fAccXYZ[0]*9.8f/1000.0f;   //單位會是 m/s^2, 每平方秒n米
    xSensorData.m_fMag[0]=fMagXYZ[0]/1000.0;
    xSensorData.m_fGyro[0]=fGyroXYZ[0]*3.141592f/180.0f;
  
  
    xSensorData.m_fAcc[1]=-fAccXYZ[1]*9.8f/1000.0f;
    xSensorData.m_fMag[1]=-fMagXYZ[2]/1000.0;
    xSensorData.m_fGyro[1]=-fGyroXYZ[1]*3.141592f/180.0f;
  
    xSensorData.m_fAcc[2]=-fAccXYZ[2]*9.8f/1000.0f;
    xSensorData.m_fMag[2]=-fMagXYZ[1]/1000.0;
    xSensorData.m_fGyro[2]=-fGyroXYZ[2]*3.141592f/180.0f;
    
  }
  
  iNEMO_AHRS_Update(&xSensorData, &xEulerAngles, &xQuat);
  */
  if(iter++ == 16)
  {
    t++;

    for(i=0; i<10; i++)
    {
    	DestnaStrLen[i]=0;
    	DestnaStrX[i]=0;
    	DestnaStrY[i]=0;
    	DestnaStrZ[i]=0;
      DestnaStrT[i]=0;
      DestnaStrStep[i]=0;
      DestnaStrLengthAvg[i]=0;
    }
    // Read data
    fAccXYZ[0] = toG((float)Axes_Data.X) * 1000; //mg單位
    fAccXYZ[1] = toG((float)Axes_Data.Y) * 1000;
    fAccXYZ[2] = toG((float)Axes_Data.Z) * 1000;   	
  	
    /* Print data via Virtual Com */
    //printf("ACC(mg): X=%f, Y=%f, Z=%f\n\r", fAccXYZ[0], fAccXYZ[1], fAccXYZ[2]);
    
    /*
    roll  = xEulerAngles.m_fRoll* 180.0f / 3.141592f;
    pitch = xEulerAngles.m_fPitch * 180.0f / 3.141592f;
    yaw   = xEulerAngles.m_fYaw * 180.0f / 3.141592f;
    */

    //printf("Attitude(deg): R=%.3f, P=%.3f, Y=%.3f\n\n\r", roll, pitch, yaw);
    iter=0;
    //============================================================//
    Length = sqrt(fAccXYZ[0]*fAccXYZ[0]+fAccXYZ[1]*fAccXYZ[1]+fAccXYZ[2]*fAccXYZ[2]);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    float2String(&DestnaStrX[0], fAccXYZ[0]);
    float2String(&DestnaStrY[0], fAccXYZ[1]);
    float2String(&DestnaStrZ[0], fAccXYZ[2]);
    float2String(&DestnaStrLen[0], Length);
    float2String(&DestnaStrT[0], t);
    float2String(&DestnaStrStep[0], step);
    float2String(&DestnaStrLengthAvg[0], LengthAvg);

    //printf("ACC(mg): X=%f, Y=%f, Z=%f Length=%.3f\n", fAccXYZ[0], fAccXYZ[1], fAccXYZ[2], Length);

    
    //=======================================================================//
    //My_Usart2_Printf("T = ");
    //printf("T = ");
    
    for(i=9; i>=0; i--)
    { 
      if(DestnaStrT[i]==0)
      {
      }	
      else
      {	
        USART_SendData(USART2, DestnaStrT[i]);
        //printf("%c", DestnaStrT[i]);
        
        while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
      }  
    }
    
    //My_Usart2_Printf("; ");
    My_Usart2_Printf(" ");
    //printf("; ");
    
    //My_Usart2_Printf("X = ");
    //printf("X = ");
    
    if(fAccXYZ[0]<0)
    {
    	USART_SendData(USART2, DestnaStrX[0]);
      //printf("%c", DestnaStrX[0]);
      
      while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
      
      for(i=6; i>0; i--)
      {
        if(DestnaStrX[i]==0)
        {
        }
        else
        {	
          USART_SendData(USART2, DestnaStrX[i]);
          //printf("%c", DestnaStrX[i]);
          
          while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
        }  
      }      
    }
    else
    {
      for(i=6; i>=0; i--)
      {
        if(DestnaStrX[i]==0)
        {
        }
        else
        {	
          USART_SendData(USART2, DestnaStrX[i]);
          //printf("%c", DestnaStrX[i]);
          
          while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
        }  
      }
    }       
    
    //My_Usart2_Printf("; ");
    My_Usart2_Printf(" ");
    //printf("; ");
    
    //My_Usart2_Printf("Y = ");
    //printf("Y = ");
    
    if(fAccXYZ[1]<0)
    {
      USART_SendData(USART2, DestnaStrY[0]);
      //printf("%c", DestnaStrY[0]);
      
      while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
      
      for(i=6; i>0; i--)
      {
        if(DestnaStrY[i]==0)
        {
        }
        else
        {	
          USART_SendData(USART2, DestnaStrY[i]);
          //printf("%c", DestnaStrY[i]);
          
          while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
        }  
      } 
    }
    else
    {	
      for(i=6; i>=0; i--)
      {
        if(DestnaStrY[i]==0)
        {
        }
        else
        {	
          USART_SendData(USART2, DestnaStrY[i]);
          //printf("%c", DestnaStrY[i]);
          
          while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
        }  
      }
    }        
    
    //My_Usart2_Printf("; ");
    My_Usart2_Printf(" ");
    //printf("; ");
    
    //My_Usart2_Printf("Z = ");
    //printf("Z = ");

    if(fAccXYZ[2]<0)
    {
      USART_SendData(USART2, DestnaStrZ[0]);
      //printf("%c", DestnaStrZ[0]);
      
      while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
      
      for(i=6; i>0; i--)
      {
        if(DestnaStrZ[i]==0)
        {
        }
        else
        {	
          USART_SendData(USART2, DestnaStrZ[i]);
          //printf("%c", DestnaStrZ[i]);
          
          while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
        }  
      }     	
    }
    else
    {	    
      for(i=6; i>=0; i--)
      {
        if(DestnaStrZ[i]==0)
        {
        }
        else
        {	
          USART_SendData(USART2, DestnaStrZ[i]);
          //printf("%c", DestnaStrZ[i]);
          
          while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
        }  
      }
    } 
    //My_Usart2_Printf("; ");
    My_Usart2_Printf(" ");
    //printf("; ");
    
    //My_Usart2_Printf("Len = ");
    //printf("Len = ");
    
    for(i=6; i>=0; i--)
    {
      if(DestnaStrLen[i]==0)
      {
      }
      else
      {	
        USART_SendData(USART2, DestnaStrLen[i]);
        //printf("%c", DestnaStrLen[i]);
        
        while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
      }  
    } 
    
    My_Usart2_Printf(" ");
    
    for(i=9; i>=0; i--)
    { 
      if(DestnaStrStep[i]==0)
      {
      }	
      else
      {	
        USART_SendData(USART2, DestnaStrStep[i]);
        //printf("%c", DestnaStrT[i]);
        
        while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
      }  
    }
    
    My_Usart2_Printf(" ");
    
    for(i=9; i>=0; i--)
    { 
      if(DestnaStrLengthAvg[i]==0)
      {
      }	
      else
      {	
        USART_SendData(USART2, DestnaStrLengthAvg[i]);
        //printf("%c", DestnaStrT[i]);
        
        while( USART_GetFlagStatus( USART2, USART_FLAG_TXE ) == RESET );
      }  
    } 
    
    
    My_Usart2_Printf("\r\n");
    //printf("\r\n");    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    if(t==1)
    {	
      LengthAvg = Length;
    }
    else if(t<100)
    {
      if(Length<(LengthAvg-600))
        LengthAvg = LengthAvg*0.9 + Length*0.1;
      else if(Length>(LengthAvg+600))
        LengthAvg = LengthAvg*0.9 + Length*0.1;
      else  	
        LengthAvg = LengthAvg*0.6 + Length*0.4;  
    }
    else
    {	
      if(Length<(LengthAvg-400))
      {
        step++;
      }
      else if(Length>(LengthAvg+400))
      {
      }

      LengthAvg = LengthAvg*0.6 + Length*0.4;
    }
    
    
    /*
      方法: 當超過threshold後 再回來threshold以下 才算一次step
    */
    /*
    if((Length>LengthMax)||(Length<LengthMin))
    {
      if(!lastLengthState) //從standard到max/min
    	{
    	}
    	else //還是在max/min
    	{
    	}	
      
      lastLengthState=1;
    }
    else
    {
      if(!lastLengthState) //還是在standard
      {
      }
      else //從max/min到standard
      {
        //step++;
      }
      lastLengthState=0;
    }
    */
    //============================================================//
    //add by ycliu at 20141223
    /*
    if(iterj<7)
    {
      xcont[iterj] = fAccXYZ[0];
      zcont[iterj] = fAccXYZ[2];
      iterj++;
    }
    else
    {
      //運算平均	  	
    	for(int i=0; i<7; i++)
    	{
    	  AVGx += xcont[i];
    	  AVGz += zcont[i];
    	}
    	
    	AVGx /= 7;
    	AVGz /= 7;
    		
    	//重新update
      iterj=0;	
      xcont[iterj] = fAccXYZ[0];
      zcont[iterj] = fAccXYZ[2];
    }
      	
      	
    //計步判斷
    if((((xcont[1]-xcont[0])>0) && ((xcont[2]-xcont[1])>0) && ((xcont[4]-xcont[3])<0)) ||
      (((xcont[1]-xcont[0])>0) && ((xcont[2]-xcont[1])>0) && ((xcont[3]-xcont[2])>0)))
    {
    
      if(((xcont[4]-xcont[3])<0) && ((xcont[5]-xcont[3])<0) && ((xcont[6]-xcont[5])<0))
      {
        if((xcont[2]<AVGx) && (xcont[3]<AVGx) && (xcont[4]<AVGx) && (fAccXYZ[2]<AVGz))
        {
          step++;
        }    	
      }  		
    }
    */
    //============================================================//  
  } 
  
  return 0;
}


/* __GSENSOR_AHRS_H */

