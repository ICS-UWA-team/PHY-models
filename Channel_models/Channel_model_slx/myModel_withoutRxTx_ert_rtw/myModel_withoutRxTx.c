/*
 * File: myModel_withoutRxTx.c
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

#include "myModel_withoutRxTx.h"
#include <math.h>
#include <stdlib.h>
#define MAX_BUFFER_SIZE                32768
#ifndef ssGetFixedStepSize
#define ssGetFixedStepSize(S)          (S).stepSize
#endif                                 /* ssGetFixedStepSize */

/* Block signals and states (default storage) */
DW rtDW;

/* External inputs (root inport signals with default storage) */
ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
ExtY rtY;

/* Real-time model */
RT_MODEL rtM_;
RT_MODEL *const rtM = &rtM_;

/* Model step function for TID0 */
void myModel_withoutRxTx_step0(void)   /* Sample time: [1.6666666666666667E-6s, 0.0s] */
{
  int32_T r;
  int32_T coefPolyphaseOffset;
  int32_T nModDFactor;
  int32_T n;
  int32_T k;
  real_T rtb_Divide;
  int32_T tmp;
  int32_T tmp_0;
  int32_T i;
  real_T rtb_MatrixMultiply2_re;
  real_T rtb_MatrixMultiply2_im;
  creal_T rtb_MatrixMultiply2;
  const creal_T *rtb_FIRRateConversion_0;

  /* Update the flag to indicate when data transfers from
   *  Sample time: [1.6666666666666667E-6s, 0.0s] to Sample time: [0.05s, 0.0s]  */
  (rtM->Timing.RateInteraction.TID0_1)++;
  if ((rtM->Timing.RateInteraction.TID0_1) > 29999) {
    rtM->Timing.RateInteraction.TID0_1 = 0;
  }

  /* S-Function (discreteNLastSamples): '<S1>/S-Function ' incorporates:
   *  Inport: '<Root>/In1'
   */
  /* S-Function block: <S1>/S-Function  */
  {
    int nSamples = 9624 ;
    int io = 0;
    int iv;
    int ind_Ps = rtDW.SFunction_IWORK.indPs;

    /* Input present value(s) */
    ((real_T *)rtDW.SFunction_PWORK.uBuffers)[ind_Ps] = rtU.In1;

    /* Output past value(s) */
    /* Output from present sample index to 0 */
    for (iv = ind_Ps; iv >= 0; --iv)
      (&rtDW.SFunction[0])[io++] = ((real_T *)rtDW.SFunction_PWORK.uBuffers)[iv];

    /* Output from end of buffer to present sample index excl. */
    for (iv = nSamples-1; iv > ind_Ps; --iv)
      (&rtDW.SFunction[0])[io++] = ((real_T *)rtDW.SFunction_PWORK.uBuffers)[iv];

    /* Update ring buffer index */
    if (++(rtDW.SFunction_IWORK.indPs) == nSamples)
      rtDW.SFunction_IWORK.indPs = 0;
  }

  /* RateTransition: '<Root>/TmpRTBAtVariable SelectorInport2' */
  if (rtM->Timing.RateInteraction.TID0_1 == 1) {
    rtDW.TmpRTBAtVariableSelectorInport2 = rtDW.TmpRTBAtVariableSelectorInpor_m;
  }

  /* End of RateTransition: '<Root>/TmpRTBAtVariable SelectorInport2' */

  /* S-Function (sdspperm2): '<Root>/Variable Selector' incorporates:
   *  Constant: '<Root>/Constant'
   */
  i = (int32_T)floor(rtDW.TmpRTBAtVariableSelectorInport2) - 1;
  if (i < 0) {
    i = 0;
  } else {
    if (i >= 3603) {
      i = 3602;
    }
  }

  memcpy(&rtDW.MatrixMultiply2[0], &rtConstP.Constant_Value[i * 401], 401U *
         sizeof(creal_T));

  /* End of S-Function (sdspperm2): '<Root>/Variable Selector' */

  /* S-Function (sdspcount2): '<Root>/Counter' */
  rtb_Divide = rtDW.Counter_Count;
  if (rtDW.Counter_Count < 30000) {
    rtDW.Counter_Count++;
  } else {
    rtDW.Counter_Count = 0U;
  }

  /* End of S-Function (sdspcount2): '<Root>/Counter' */

  /* RateTransition: '<Root>/TmpRTBAtVariable Selector1Inport2' */
  if (rtM->Timing.RateInteraction.TID0_1 == 1) {
    rtDW.TmpRTBAtVariableSelector1Inport = rtDW.TmpRTBAtVariableSelector1Inpo_c;
  }

  /* End of RateTransition: '<Root>/TmpRTBAtVariable Selector1Inport2' */

  /* S-Function (sdspperm2): '<Root>/Variable Selector1' */
  i = (int32_T)floor(rtDW.TmpRTBAtVariableSelector1Inport) - 1;
  if (i < 0) {
    i = 0;
  } else {
    if (i >= 3603) {
      i = 3602;
    }
  }

  /* Product: '<Root>/Divide' incorporates:
   *  Constant: '<Root>/Constant2'
   */
  rtb_Divide /= 30000.0;
  for (r = 0; r < 401; r++) {
    /* S-Function (sdspperm2): '<Root>/Variable Selector1' incorporates:
     *  Constant: '<Root>/Constant'
     */
    rtb_FIRRateConversion_0 = &rtConstP.Constant_Value[401 * i + r];

    /* Product: '<Root>/Matrix Multiply1' incorporates:
     *  Constant: '<Root>/Constant1'
     *  Sum: '<Root>/Minus'
     */
    rtb_MatrixMultiply2_re = (1.0 - rtb_Divide) * rtDW.MatrixMultiply2[r].re;
    rtb_MatrixMultiply2_im = (1.0 - rtb_Divide) * rtDW.MatrixMultiply2[r].im;
    rtDW.MatrixMultiply[r].re = rtb_MatrixMultiply2_re;
    rtDW.MatrixMultiply[r].im = rtb_MatrixMultiply2_im;

    /* Sum: '<Root>/Minus1' incorporates:
     *  Product: '<Root>/Matrix Multiply2'
     */
    rtDW.MatrixMultiply[r].re = rtb_FIRRateConversion_0->re * rtb_Divide +
      rtb_MatrixMultiply2_re;
    rtDW.MatrixMultiply[r].im = rtb_FIRRateConversion_0->im * rtb_Divide +
      rtb_MatrixMultiply2_im;

    /* Product: '<Root>/Matrix Multiply1' */
    rtb_MatrixMultiply2.re = rtb_MatrixMultiply2_re;
    rtb_MatrixMultiply2.im = rtb_MatrixMultiply2_im;

    /* S-Function (sdspperm2): '<Root>/Variable Selector1' */
    rtDW.FIRRateConversion[r] = *rtb_FIRRateConversion_0;

    /* Product: '<Root>/Matrix Multiply1' */
    rtDW.MatrixMultiply2[r] = rtb_MatrixMultiply2;
  }

  /* S-Function (sdspupfirdn2): '<Root>/FIR Rate Conversion' */
  i = 0;
  r = 0;

  /* Update inBufIdx and inputChannelOffset for current channel */
  tmp = rtDW.FIRRateConversion_InBufIdx;
  for (n = 0; n < 401; n++) {
    nModDFactor = n % 401;

    /* Read input into inBufArray */
    rtDW.FIRRateConversion_InBuf[tmp] = rtDW.MatrixMultiply[r];
    r++;

    /* Generate outputs (if any) for current input n */
    for (k = rtConstP.FIRRateConversion_StartIdx[nModDFactor]; k <
         rtConstP.FIRRateConversion_StopIdx[nModDFactor]; k++) {
      rtb_Divide = 0.0;
      rtb_MatrixMultiply2_re = 0.0;
      coefPolyphaseOffset = rtConstP.FIRRateConversion_PolyphaseSele[k] * 24;
      for (tmp_0 = tmp; tmp_0 < 24; tmp_0++) {
        rtb_MatrixMultiply2_im = rtConstP.FIRRateConversion_FILTER
          [(coefPolyphaseOffset + tmp_0) - tmp];
        rtb_Divide += rtb_MatrixMultiply2_im *
          rtDW.FIRRateConversion_InBuf[tmp_0].re;
        rtb_MatrixMultiply2_re += rtb_MatrixMultiply2_im *
          rtDW.FIRRateConversion_InBuf[tmp_0].im;
      }

      for (tmp_0 = 0; tmp_0 < tmp; tmp_0++) {
        rtb_MatrixMultiply2_im = rtConstP.FIRRateConversion_FILTER
          [((coefPolyphaseOffset - tmp) + tmp_0) + 24];
        rtb_Divide += rtb_MatrixMultiply2_im *
          rtDW.FIRRateConversion_InBuf[tmp_0].re;
        rtb_MatrixMultiply2_re += rtb_MatrixMultiply2_im *
          rtDW.FIRRateConversion_InBuf[tmp_0].im;
      }

      rtDW.FIRRateConversion[i].re = rtb_Divide;
      rtDW.FIRRateConversion[i].im = rtb_MatrixMultiply2_re;
      i++;
    }

    /* Decrement inBufIdx, wrap if necessary */
    if (tmp == 0) {
      tmp = 23;
    } else {
      tmp--;
    }
  }

  /* Update inBufIdx */
  rtDW.FIRRateConversion_InBufIdx = tmp;

  /* End of S-Function (sdspupfirdn2): '<Root>/FIR Rate Conversion' */

  /* Sum: '<Root>/Sum of Elements' */
  rtb_Divide = -0.0;
  rtb_MatrixMultiply2_re = 0.0;
  for (i = 0; i < 9624; i++) {
    /* Sum: '<Root>/Sum of Elements' incorporates:
     *  Product: '<Root>/Matrix Multiply'
     */
    rtb_Divide += rtDW.SFunction[i] * rtDW.FIRRateConversion[i].re;
    rtb_MatrixMultiply2_re += rtDW.SFunction[i] * rtDW.FIRRateConversion[i].im;
  }

  /* Outport: '<Root>/Out1' incorporates:
   *  Sum: '<Root>/Sum of Elements'
   */
  rtY.Out1.re = rtb_Divide;
  rtY.Out1.im = rtb_MatrixMultiply2_re;
}

/* Model step function for TID1 */
void myModel_withoutRxTx_step1(void)   /* Sample time: [0.05s, 0.0s] */
{
  int32_T rtb_Minus2;
  int32_T rtb_Minus3;

  /* S-Function (sdspcount2): '<Root>/Counter1' */
  rtb_Minus2 = rtDW.Counter1_Count;
  if (rtDW.Counter1_Count < 3602) {
    rtDW.Counter1_Count++;
  } else {
    rtDW.Counter1_Count = 0U;
  }

  /* End of S-Function (sdspcount2): '<Root>/Counter1' */

  /* S-Function (sdspcount2): '<Root>/Counter2' */
  rtb_Minus3 = rtDW.Counter2_Count;
  if (rtDW.Counter2_Count < 3602) {
    rtDW.Counter2_Count++;
  } else {
    rtDW.Counter2_Count = 0U;
  }

  /* End of S-Function (sdspcount2): '<Root>/Counter2' */

  /* Sum: '<Root>/Minus3' incorporates:
   *  Constant: '<Root>/Constant1'
   */
  rtb_Minus3++;

  /* Sum: '<Root>/Minus2' incorporates:
   *  Constant: '<Root>/Constant1'
   */
  rtb_Minus2++;

  /* Update for RateTransition: '<Root>/TmpRTBAtVariable SelectorInport2' */
  rtDW.TmpRTBAtVariableSelectorInpor_m = rtb_Minus2;

  /* Update for RateTransition: '<Root>/TmpRTBAtVariable Selector1Inport2' */
  rtDW.TmpRTBAtVariableSelector1Inpo_c = rtb_Minus3;
}

/* Model initialize function */
void myModel_withoutRxTx_initialize(void)
{
  /* Start for S-Function (discreteNLastSamples): '<S1>/S-Function ' incorporates:
   *  Inport: '<Root>/In1'
   */

  /* S-Function Block: <S1>/S-Function  */
  {
    static real_T dnls_buffer[1 * 9624];
    rtDW.SFunction_PWORK.uBuffers = (void *)&dnls_buffer[0];
  }

  {
    int ids;

    /* Assign default sample(s) */
    if (rtDW.SFunction_PWORK.uBuffers != NULL) {
      for (ids = 0; ids < 9624; ++ids)
        ((real_T *)rtDW.SFunction_PWORK.uBuffers)[ids] = (real_T)0.0;
    }

    /* Set work values */
    rtDW.SFunction_IWORK.indPs = 0;
  }

  /* InitializeConditions for S-Function (sdspcount2): '<Root>/Counter2' */
  rtDW.Counter2_Count = 1U;
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
