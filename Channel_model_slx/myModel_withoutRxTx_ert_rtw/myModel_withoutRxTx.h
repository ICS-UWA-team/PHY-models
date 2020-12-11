/*
 * File: myModel_withoutRxTx.h
 *
 * Code generated for Simulink model 'myModel_withoutRxTx'.
 *
 * Model version                  : 1.101
 * Simulink Coder version         : 9.0 (R2018b) 24-May-2018
 * C/C++ source code generated on : Fri Nov 13 17:06:26 2020
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM Cortex
 * Emulation hardware selection:
 *    Differs from embedded hardware (Custom Processor->MATLAB Host Computer)
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_myModel_withoutRxTx_h_
#define RTW_HEADER_myModel_withoutRxTx_h_
#include <math.h>
#include <string.h>
#ifndef myModel_withoutRxTx_COMMON_INCLUDES_
# define myModel_withoutRxTx_COMMON_INCLUDES_
#include "rtwtypes.h"
#endif                                 /* myModel_withoutRxTx_COMMON_INCLUDES_ */

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

/* Forward declaration for rtModel */
typedef struct tag_RTM RT_MODEL;

/* Block signals and states (default storage) for system '<Root>' */
typedef struct {
  creal_T FIRRateConversion_InBuf[24]; /* '<Root>/FIR Rate Conversion' */
  creal_T MatrixMultiply[9624];        /* '<Root>/Matrix Multiply' */
  creal_T FIRRateConversion[9624];     /* '<Root>/FIR Rate Conversion' */
  creal_T MatrixMultiply2[401];        /* '<Root>/Matrix Multiply2' */
  real_T SFunction[9624];              /* '<S1>/S-Function ' */
  real_T TmpRTBAtVariableSelectorInport2;/* '<Root>/Minus2' */
  real_T TmpRTBAtVariableSelector1Inport;/* '<Root>/Minus3' */
  real_T TmpRTBAtVariableSelectorInpor_m;/* synthesized block */
  real_T TmpRTBAtVariableSelector1Inpo_c;/* synthesized block */
  struct {
    void *uBuffers;
  } SFunction_PWORK;                   /* '<S1>/S-Function ' */

  struct {
    int_T indPs;
    int_T bufSz;
  } SFunction_IWORK;                   /* '<S1>/S-Function ' */

  int32_T FIRRateConversion_InBufIdx;  /* '<Root>/FIR Rate Conversion' */
  uint16_T Counter_Count;              /* '<Root>/Counter' */
  uint16_T Counter1_Count;             /* '<Root>/Counter1' */
  uint16_T Counter2_Count;             /* '<Root>/Counter2' */
} DW;

/* Constant parameters (default storage) */
typedef struct {
  /* Expression: a.h
   * Referenced by: '<Root>/FIR Rate Conversion'
   */
  real_T FIRRateConversion_FILTER[230976];

  /* Expression: hmat
   * Referenced by: '<Root>/Constant'
   */
  creal_T Constant_Value[1444803];

  /* Computed Parameter: FIRRateConversion_PolyphaseSele
   * Referenced by: '<Root>/FIR Rate Conversion'
   */
  int32_T FIRRateConversion_PolyphaseSele[9624];

  /* Computed Parameter: FIRRateConversion_StartIdx
   * Referenced by: '<Root>/FIR Rate Conversion'
   */
  int32_T FIRRateConversion_StartIdx[401];

  /* Computed Parameter: FIRRateConversion_StopIdx
   * Referenced by: '<Root>/FIR Rate Conversion'
   */
  int32_T FIRRateConversion_StopIdx[401];
} ConstP;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T In1;                          /* '<Root>/In1' */
} ExtU;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  creal_T Out1;                        /* '<Root>/Out1' */
} ExtY;

/* Real-time Model Data Structure */
struct tag_RTM {
  const char_T * volatile errorStatus;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    struct {
      uint16_T TID0_1;
    } RateInteraction;
  } Timing;
};

/* Block signals and states (default storage) */
extern DW rtDW;

/* External inputs (root inport signals with default storage) */
extern ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY rtY;

/* Constant parameters (default storage) */
extern const ConstP rtConstP;

/* Model entry point functions */
extern void myModel_withoutRxTx_initialize(void);
extern void myModel_withoutRxTx_step0(void);
extern void myModel_withoutRxTx_step1(void);

/* Real-time Model object */
extern RT_MODEL *const rtM;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'myModel_withoutRxTx'
 * '<S1>'   : 'myModel_withoutRxTx/Discrete Shift Register'
 */
#endif                                 /* RTW_HEADER_myModel_withoutRxTx_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
