/*	____________________________________________________________________
 *
 *	This is an example sketch for Linux OS, that shows how to control
 *	SimpleBGC-driven gimbal via Serial API. API specs are available at
 *	http://www.basecamelectronics.com/serialapi/
 *	____________________________________________________________________
 */

#include "sbgc32.h"

#define TIMESTEP_WAIT 10
#define TIMESTEP_INC 0.0002
#define LIMIT_LOW -45
#define LIMIT_HIGH 0

/* ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾ */
/*   					Global Software Objects  					  */
/* __________________________________________________________________ */

GeneralSBGC_t SBGC32_Device;

static Control_t Control;
static ControlConfig_t ControlConfig;

static BoardInfo_t BoardInfo;
static BoardInfo3_t BoardInfo3;
static MainParams3_t MainParams3;
static MainParamsExt_t MainParamsExt;
static MainParamsExt2_t MainParamsExt2;
static MainParamsExt3_t MainParamsExt3;

static RealTimeDataCustom_t RealTimeDataCustom;
static RealTimeData_t RealTimeData;

static AdjVarGeneral_t AdjVarGeneral[3];

static DataStreamInterval_t DataStreamInterval;

static BeeperSettings_t BeeperSettings;

// static      ConfirmationState_t Confirm;

static ui8 DataStreamBuff[20];

TxRxStatus_t PrintBoardParameters(Profile_t slot);
TxRxStatus_t SBGC32_DemoControl(void);
void PrintDataStream(ui8 *pBuff);

/*  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = */

typedef struct
{
    double Y;
    double Yd;
    double Ydd;
    double sampledTime;
} Trap_Step_t;

typedef struct
{
    double vel_limit;   // [count/s]
    double accel_limit; // [count/s^2]
    double decel_limit; // [count/s^2]
                        // float A_per_css  ;      // [A/(count/s^2)]
    double Xi_;
    double Xf_;
    double Vi_;

    double Ar_;
    double Vr_;
    double Dr_;

    double Ta_;
    double Tv_;
    double Td_;
    double Tf_;

    double yAccel_;
    double lastTimeStep;
    uint8_t shouldExecute;
} Trap_Config_t;

#define SQ(x) ((x) * (x))
uint8_t planTrapezoidal(Trap_Config_t *Cfg, double Xf, double Xi, float Vi,
                        float Vmax, float Amax, float Dmax);

uint8_t evalTrapezoidal(Trap_Config_t *Cfg, double t, Trap_Step_t *Step_t);
float sign_hard(float val)
{
    return (signbit(val)) ? -1.0f : 1.0f;
}
volatile Trap_Step_t TStep[3];
volatile Trap_Config_t TCfg[3];
volatile double v_pos[6];

double dmap(double input, double input_start, double input_end, double output_start, double output_end)
{

    double slope;

    if ((input < input_start) || (input > input_end))
        return 0;

    slope = 1.0 * (output_end - output_start) / (input_end - input_start);

    // return  output_start + ddround(slope * (input - input_start));
    return output_start + (slope * (input - input_start));
}

double clamp(double x, double lowerlimit, double upperlimit)
{
    if (x < lowerlimit)
        x = lowerlimit;
    if (x > upperlimit)
        x = upperlimit;
    return x;
}

double smootherstep(double edge0, double edge1, double x)
{
    // Scale, and clamp x to 0..1 range
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    // Evaluate polynomial
    return x * x * x * (x * (x * 6 - 15) + 10);
    // return  ((x) * (x) * (3 - 2 * (x)));
}

double getIntValue(double start, double finish, double time, double totalTime)
{

    double t;
    double sp;
    double ep;
    double ret;
    t = dmap(time, 0, totalTime, 0, 1);

    ret = smootherstep(0, 1, t);
    return dmap(ret, 0, 1, start, finish); // ret = map(ret,)
}

// Symbol                     Description
// Ta, Tv and Td              Duration of the stages of the AL profile
// Xi and Vi                  Adapted initial conditions for the AL profile
// Xf                         Position set-point
// s                          Direction (sign) of the trajectory
// Vmax, Amax, Dmax and jmax  Kinematic bounds
// Ar, Dr and Vr              Reached values of acceleration and velocity

uint8_t planTrapezoidal(Trap_Config_t *Cfg, double Xf, double Xi, float Vi,
                        float Vmax, float Amax, float Dmax)
{

    double dX = Xf - Xi;                          // Distance to travel
    double stop_dist = (Vi * Vi) / (2.0f * Dmax); // Minimum stopping distance
    float dXstop = copysignf(stop_dist, Vi);      // Minimum stopping displacement
    float s = sign_hard(dX - dXstop);             // Sign of coast velocity (if any)
    Cfg->Ar_ = s * Amax;                          // Maximum Acceleration (signed)
    Cfg->Dr_ = -s * Dmax;                         // Maximum Deceleration (signed)
    Cfg->Vr_ = s * Vmax;                          // Maximum Velocity (signed)

    // If we start with a speed faster than cruising, then we need to decel instead of accel
    // aka "double deceleration move" in the paper
    if ((s * Vi) > (s * Cfg->Vr_))
    {
        Cfg->Ar_ = -s * Amax;
    }

    // Time to accel/decel to/from Vr (cruise speed)
    Cfg->Ta_ = (Cfg->Vr_ - Vi) / Cfg->Ar_;
    Cfg->Td_ = -Cfg->Vr_ / Cfg->Dr_;

    // Integral of velocity ramps over the full accel and decel times to get
    // minimum displacement required to reach cuising speed
    float dXmin = 0.5f * Cfg->Ta_ * (Cfg->Vr_ + Vi) + 0.5f * Cfg->Td_ * Cfg->Vr_;

    // Are we displacing enough to reach cruising speed?
    if (s * dX < s * dXmin)
    {
        // Short move (triangle profile)
        Cfg->Vr_ = s * sqrtf(fmax((Cfg->Dr_ * SQ(Vi) + 2 * Cfg->Ar_ * Cfg->Dr_ * dX) / (Cfg->Dr_ - Cfg->Ar_), 0.0f));
        Cfg->Ta_ = fmax(0.0f, (Cfg->Vr_ - Vi) / Cfg->Ar_);
        Cfg->Td_ = fmax(0.0f, -Cfg->Vr_ / Cfg->Dr_);
        Cfg->Tv_ = 0.0f;
    }
    else
    {
        // Long move (trapezoidal profile)
        Cfg->Tv_ = (dX - dXmin) / Cfg->Vr_;
    }

    // Fill in the rest of the values used at evaluation-time
    Cfg->Tf_ = Cfg->Ta_ + Cfg->Tv_ + Cfg->Td_;
    Cfg->Xi_ = Xi;
    Cfg->Xf_ = Xf;
    Cfg->Vi_ = Vi;
    Cfg->yAccel_ = Xi + Vi * Cfg->Ta_ + 0.5f * Cfg->Ar_ * SQ(Cfg->Ta_); // pos at end of accel phase

    Cfg->lastTimeStep = 0;
    Cfg->shouldExecute = 1;
    return 1;
}

uint8_t evalTrapezoidal(Trap_Config_t *Cfg, double t, Trap_Step_t *Step_t)
{
    if (t < 0.0f)
    { // Initial Condition
        Step_t->Y = Cfg->Xi_;
        Step_t->Yd = Cfg->Vi_;
        Step_t->Ydd = 0.0f;
    }
    else if (t < Cfg->Ta_)
    { // Accelerating
        Step_t->Y = Cfg->Xi_ + Cfg->Vi_ * t + 0.5f * Cfg->Ar_ * SQ(t);
        Step_t->Yd = Cfg->Vi_ + Cfg->Ar_ * t;
        Step_t->Ydd = Cfg->Ar_;
    }
    else if (t < Cfg->Ta_ + Cfg->Tv_)
    { // Coasting
        Step_t->Y = Cfg->yAccel_ + Cfg->Vr_ * (t - Cfg->Ta_);
        Step_t->Yd = Cfg->Vr_;
        Step_t->Ydd = 0.0f;
    }
    else if (t < Cfg->Tf_)
    { // Deceleration
        float td = t - Cfg->Tf_;
        Step_t->Y = Cfg->Xf_ + 0.5f * Cfg->Dr_ * SQ(td);
        Step_t->Yd = Cfg->Dr_ * td;
        Step_t->Ydd = Cfg->Dr_;
    }
    else if (t >= Cfg->Tf_)
    { // Final Condition
        Step_t->Y = Cfg->Xf_;
        Step_t->Yd = 0.0f;
        Step_t->Ydd = 0.0f;
    }
    else
    {
        // TODO: report error here
    }
    Step_t->sampledTime = t;
    return 0;
}

void init_gimbal()
{
    // SBGC_1.Drv = malloc(sizeof(Driver_t));
    // std::cout << "Connecting to serial port: " << SBGC_SERIAL_PORT<< "\n";
    // DriverInit(SBGC_1.Drv, SBGC_SERIAL_PORT);

    /* High Layer Init */
    // SBGC32_DefaultInit(&SBGC_1, PortTransmitData, PortReceiveByte, GetAvailableBytes, PrintDebugData, GetTimeMs, SBGC_PROTOCOL_V2);

    /* - - - - - - - - - High Layer Software Init - - - - - - - - - - */

    // ControlConfig.AxisCC[ROLL].angleLPF  = 6;
    // ControlConfig.AxisCC[PITCH].angleLPF = 7;
    ControlConfig.AxisCC[YAW].angleLPF = 7;

    // ControlConfig.AxisCC[ROLL].angleLPF  = 6;
    // ControlConfig.AxisCC[PITCH].speedLPF = 7;
    ControlConfig.AxisCC[YAW].speedLPF = 7;
    ControlConfig.flags = RTCCF_CONTROL_CONFIG_FLAG_NO_CONFIRM;

    // Control.controlMode[ROLL]  = CtrlM_MODE_ANGLE | CtrlF_CONTROL_FLAG_TARGET_PRECISE |CtrlF_CONTROL_FLAG_HIGH_RES_SPEED;
    // Control.controlMode[PITCH] = CtrlM_MODE_ANGLE | CtrlF_CONTROL_FLAG_TARGET_PRECISE;
    Control.controlMode[YAW] = CtrlM_MODE_ANGLE | CtrlF_CONTROL_FLAG_TARGET_PRECISE | CtrlF_CONTROL_FLAG_HIGH_RES_SPEED;

    // Control.AxisC[ROLL].angle  = 0;
    // Control.AxisC[PITCH].angle = 0;
    Control.AxisC[YAW].angle = 0;

    // Control.AxisC[PITCH].speed = 500;
    Control.AxisC[YAW].speed = 500;

    /* SBGC32_Reset(&SBGC_1, RF_RESET_WITH_RESTORING_STATES, 500);
     *     SBGC32_CheckConfirmation(&SBGC_1, &Confirm, CMD_RESET);
     *         sleep(5); */

    //  SBGC32_ControlConfig(&SBGC32_Device, &ControlConfig, &Confirm);
    SBGC32_ControlConfig(&SBGC32_Device, &ControlConfig);

    PrintBoardParameters(P_CURRENT_PROFILE);
}

// LPF 0..15
// accel 0..65535 degs/sec^2

// control modes
// CtrlM_MODE_NO_CONTROL
// CtrlM_MODE_SPEED
// CtrlM_MODE_ANGLE
// CtrlM_MODE_SPEED_ANGLE
// CtrlM_MODE_RC
// CtrlM_MODE_ANGLE_REL_FRAME
// CtrlM_MODE_RC_HIGH_RES

// flags
// CtrlF_CONTROL_FLAG_TARGET_PRECISE
// CtrlF_CONTROL_FLAG_AUTO_TASK
// CtrlF_CONTROL_FLAG_FORCE_RC_SPEED
// CtrlF_CONTROL_FLAG_HIGH_RES_SPEED

void set_gimbal_pitch_yaw(double pitch_deg, double yaw_deg)
{
    Control.AxisC[YAW].angle = DEGREE_TO_ANGLE_INT(yaw_deg);
    Control.AxisC[PITCH].angle = DEGREE_TO_ANGLE_INT(pitch_deg);
    SBGC32_Control(&SBGC32_Device, &Control);
}
void set_gimbal_yaw(double yaw_deg)
{
    Control.AxisC[YAW].angle = DEGREE_TO_ANGLE_INT(yaw_deg);
    SBGC32_Control(&SBGC32_Device, &Control);
}
void set_gimbal_yaw_PSA(double deg, double speed, double accel)
{
    Control.AxisC[YAW].angle = DEGREE_TO_ANGLE_INT(deg);
    Control.AxisC[YAW].speed = speed; // CtrlF_CONTROL_FLAG_HIGH_RES_SPEED 0.001 deg/sec or 0.1220740379 deg/sec
    // ControlConfig.AxisCC[YAW].AccLimit = accel;
    SBGC32_Control(&SBGC32_Device, &Control);
    // SBGC32_ControlConfig(&SBGC32_Device, &ControlConfig);
}

double ddround(double d)
{
    return floor(d + 0.5);
}

double map(double input, double input_start, double input_end, double output_start, double output_end)
{

    double slope;

    if ((input < input_start) || (input > input_end))
        return 0;

    slope = 1.0 * (output_end - output_start) / (input_end - input_start);

    // return  output_start + ddround(slope * (input - input_start));
    return output_start + (slope * (input - input_start));
}

void run_sinuosoidal(void)
{

    double timeStep = 0.0f;
    double yaw_pos = 0.0f;

    while (1)
    {
        yaw_pos = map(sin(timeStep), -1, +1, LIMIT_LOW, LIMIT_HIGH);
        set_gimbal_yaw(yaw_pos);
        timeStep = timeStep + TIMESTEP_INC;
        DELAY_MS_(TIMESTEP_WAIT - 1);
        if (fmod(timeStep, 2.0f))
        {
            printf("requested position = %f\n", yaw_pos);
        }
    }
}

void run_trap(void)
{
    unsigned char breakWhile = 1;
    unsigned char CHAN = 0;
    double timeStep = 0.0f;
    double yaw_pos = 0.0f;

    double start_position = LIMIT_LOW;
    double required_position = LIMIT_HIGH;
    double start_velocity = 0;
    double max_velocity = 20.0f; // 0.001 degs /sec -> max is 65535 -> 65.5 degs sec
    double max_acceleration = 20.0f;

    // move to init and verify the position or just wait

    printf("Moving to home position.\n");
    set_gimbal_yaw_PSA(0, 10, 0);
    DELAY_MS_(8000);
    printf("SmoothStep(fx) Trapezoidal started.\n");

    planTrapezoidal((Trap_Config_t *)&TCfg[CHAN], required_position, start_position,
                    start_velocity, max_velocity, max_acceleration,
                    max_acceleration);

    while (breakWhile)
    {

        evalTrapezoidal((Trap_Config_t *)&TCfg[CHAN], (float)((float)TCfg[CHAN].lastTimeStep / 1000.0f),
                        (Trap_Step_t *)&TStep[CHAN]);
        if ((TCfg[CHAN].lastTimeStep / 1000) < TCfg[CHAN].Ta_)
        {
            TStep[CHAN].Yd = getIntValue(0, TCfg[CHAN].Vr_, TCfg[CHAN].lastTimeStep - 1, (TCfg[CHAN].Ta_ * 1000) - 1);
            if (TCfg[CHAN].lastTimeStep == 0)
            {
                v_pos[CHAN] = TCfg[CHAN].Xi_;
            }
            else
            {
                v_pos[CHAN] = v_pos[CHAN] + (TStep[CHAN].Yd / 1000.0f);
            }
        }
        else if ((TCfg[CHAN].lastTimeStep / 1000) < TCfg[CHAN].Ta_ + TCfg[CHAN].Tv_)
        {
            // Coasting
            v_pos[CHAN] = TStep[CHAN].Y;
        }
        else if ((TCfg[CHAN].lastTimeStep / 1000) < TCfg[CHAN].Tf_)
        {
            TStep[CHAN].Yd = getIntValue(TCfg[CHAN].Vr_, 0, (TCfg[CHAN].lastTimeStep - 1 - (TCfg[CHAN].Ta_ * 1000)) - (TCfg[CHAN].Tv_ * 1000), (TCfg[CHAN].Ta_ * 1000) - 1);
            v_pos[CHAN] = v_pos[CHAN] + (TStep[CHAN].Yd / 1000.0f);
            if (fabs(TStep[CHAN].Yd) < 1)
            {
                TStep[CHAN].Yd = 1;
            }
        }
        // set_gimbal_yaw(TStep[CHAN].Y);
        // if (TCfg[CHAN].lastTimeStep % 10==0){
        if (!fmod(TCfg[CHAN].lastTimeStep, 10))
        {
            set_gimbal_yaw_PSA(TStep[CHAN].Y, TStep[CHAN].Yd * 1000, 0);
            printf("Sent %lf %lf\n", TStep[CHAN].Y, TStep[CHAN].Yd * 1000);
        }
        if ((TStep[CHAN].Y == TCfg[CHAN].Xf_) || ((TCfg[CHAN].lastTimeStep / 1000) >= TCfg[CHAN].Tf_))
        {
            breakWhile = 0;
            printf("Finished!\n");
        }
        TCfg[CHAN].lastTimeStep++;

        DELAY_MS_(1);

    } // while
}

int main()
{
    /* ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾ */
    /*                         Initialization                         */
    /* ______________________________________________________________ */

    /* SimpleBGC32 Init */
    SBGC32_Init(&SBGC32_Device);
    init_gimbal();

    run_sinuosoidal();

#ifdef original
    /*  - - - - - - - - - - Software Initialization - - - - - - - - - */

    /* SimpleBGC32 Init */
    SBGC32_Init(&SBGC32_Device);

    /* Control Configurations */
    ControlConfig.AxisCC[ROLL].angleLPF = 6;
    ControlConfig.AxisCC[PITCH].angleLPF = 6;
    ControlConfig.AxisCC[YAW].angleLPF = 7;

    ControlConfig.AxisCC[ROLL].angleLPF = 6;
    ControlConfig.AxisCC[PITCH].speedLPF = 6;
    ControlConfig.AxisCC[YAW].speedLPF = 7;
    ControlConfig.flags = RTCCF_CONTROL_CONFIG_FLAG_NO_CONFIRM;

    Control.controlMode[ROLL] = CtrlM_MODE_ANGLE | CtrlF_CONTROL_FLAG_TARGET_PRECISE;
    Control.controlMode[PITCH] = CtrlM_MODE_ANGLE | CtrlF_CONTROL_FLAG_TARGET_PRECISE;
    Control.controlMode[YAW] = CtrlM_MODE_ANGLE | CtrlF_CONTROL_FLAG_TARGET_PRECISE;

    Control.AxisC[ROLL].angle = 0;
    Control.AxisC[PITCH].angle = 0;
    Control.AxisC[YAW].angle = 0;

    Control.AxisC[PITCH].speed = SPEED_TO_VALUE(50);
    Control.AxisC[YAW].speed = SPEED_TO_VALUE(50);

    /* Data Stream Configurations */
    DataStreamInterval.cmdID = CMD_REALTIME_DATA_CUSTOM;
    DataStreamInterval.intervalMs = 1000;
    DataStreamInterval.syncToData = STD_SYNC_OFF;

    /* For more information see the SBGC32_RequestRealTimeDataCustom function.
       Total packets length = 20 bytes:
       ui16 timestampMs						 i16 [3]				i16 [3]			i16 [3] */
    ui32 DataStreamIntervalConfig = RTDCF_STATOR_ROTOR_ANGLE | RTDCF_GYRO_DATA | RTDCF_ACC_DATA;
    memcpy(DataStreamInterval.config, &DataStreamIntervalConfig, sizeof(DataStreamIntervalConfig));

    /* Adj Vars Setting. SBGC_ADJ_VARS_REF_INFO parameter must be SET_ON  */
    InitAdjVar(&AdjVarGeneral[0], ADJ_VAL_ACC_LIMITER_ROLL);
    InitAdjVar(&AdjVarGeneral[1], ADJ_VAL_ACC_LIMITER_PITCH);
    InitAdjVar(&AdjVarGeneral[2], ADJ_VAL_ACC_LIMITER_YAW);

    /* - - - - - - - - - - - - Program Launch - - - - - - - - - - - - */

    /* SBGC32_Reset(&SBGC32_Device, RF_RESTART_CONFIRMATION, 5000);
    SBGC32_CheckConfirmation(&SBGC32_Device, CMD_RESET);
    DELAY_MS_(5000); */

    PrintBoardParameters(P_CURRENT_PROFILE);

    SBGC32_ControlConfig(&SBGC32_Device, &ControlConfig);
    SBGC32_DemoControl();

    SBGC32_RequestDataStream(&SBGC32_Device, &DataStreamInterval);

    /*  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = */

    while (1)
    {
        /* ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾ */
        /* 						  Start Worker Cycle					  */
        /* ______________________________________________________________ */

        SBGC32_ParseDataStream(&SBGC32_Device, DataStreamBuff, (SBGC_Command_t)DataStreamInterval.cmdID);
        PrintDataStream(DataStreamBuff);

        DELAY_MS_(DataStreamInterval.intervalMs - 1);

        /*  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
    }
#endif
    return 0;
}

TxRxStatus_t PrintBoardParameters(Profile_t slot)
{
    SBGC32_ReadBoardInfo(&SBGC32_Device, &BoardInfo, 0);
    SBGC32_ReadBoardInfo3(&SBGC32_Device, &BoardInfo3);

    SBGC32_ReadParams3(&SBGC32_Device, &MainParams3, slot);
    SBGC32_ReadParamsExt(&SBGC32_Device, &MainParamsExt, slot);
    SBGC32_ReadParamsExt2(&SBGC32_Device, &MainParamsExt2, slot);
    SBGC32_ReadParamsExt3(&SBGC32_Device, &MainParamsExt3, slot);

    SBGC32_ReadRealTimeData4(&SBGC32_Device, &RealTimeData);

    char boardVersionStr[4];
    char firmwareVersionStr[7];

    FormatBoardVersion(&SBGC32_Device, BoardInfo.boardVer, boardVersionStr);
    FormatFirmwareVersion(&SBGC32_Device, BoardInfo.firmwareVer, firmwareVersionStr);

    PrintMessage(&SBGC32_Device, TEXT_SIZE_((char *)"Board Version: "));
    PrintMessage(&SBGC32_Device, TEXT_SIZE_(boardVersionStr));
    PrintMessage(&SBGC32_Device, TEXT_SIZE_((char *)" \n"));
    PrintMessage(&SBGC32_Device, TEXT_SIZE_((char *)"Firmware Version: "));
    PrintMessage(&SBGC32_Device, TEXT_SIZE_(firmwareVersionStr));
    PrintMessage(&SBGC32_Device, TEXT_SIZE_((char *)" \n"));

    PrintStructElement(&SBGC32_Device, &BoardInfo3.flashSize, "Flash Size =", _UNSIGNED_CHAR_);

    PrintStructElement(&SBGC32_Device, &MainParams3.profileID + 1, "Current profile #", _UNSIGNED_CHAR_); // Note: 1 --> 5
    PrintStructElement(&SBGC32_Device, &MainParams3.AxisCMP3[ROLL].p, "Roll P =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.AxisCMP3[ROLL].i, "Roll I =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.AxisCMP3[ROLL].d, "Roll D =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.AxisCMP3[PITCH].p, "Pitch P =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.AxisCMP3[PITCH].i, "Pitch I =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.AxisCMP3[PITCH].d, "Pitch D =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.AxisCMP3[YAW].p, "Yaw P =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.AxisCMP3[YAW].i, "Yaw I =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.AxisCMP3[YAW].d, "Yaw D =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.AccLimiterAll, "Acc Limiter All = ", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.AxisRC_MP3[ROLL].RC_MaxAngle, "RC Max Angle =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &MainParams3.AxisRC_MP3[YAW].RC_MinAngle, "RC Min Angle =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &MainParams3.RC_MapROLL, "RC Map Roll =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.RC_MapPITCH, "RC Map Pitch =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.RC_MapYAW, "RC Map Yaw =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.RC_MapCmd, "RC Map Cmd =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.RC_MapFC_ROLL, "RC Map FC Roll =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParams3.RC_MapFC_PITCH, "RC Map FC Pitch =", _UNSIGNED_CHAR_);

    PrintStructElement(&SBGC32_Device, &MainParamsExt.LPF_Freq[ROLL], "LPF Frequency Roll =", _UNSIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &MainParamsExt.LPF_Freq[PITCH], "LPF Frequency Pitch =", _UNSIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &MainParamsExt.LPF_Freq[YAW], "LPF Frequency Yaw =", _UNSIGNED_SHORT_);

    PrintStructElement(&SBGC32_Device, &MainParamsExt2.frameIMU_LPF_Freq, "Frame IMU LPF Freq =", _UNSIGNED_CHAR_);
    PrintStructElement(&SBGC32_Device, &MainParamsExt2.timelapseTime, "Timelapse Time =", _UNSIGNED_SHORT_);

    PrintStructElement(&SBGC32_Device, &MainParamsExt3.motorStartupDelay, "Motor Startup Delay =", _UNSIGNED_SHORT_);

    PrintMessage(&SBGC32_Device, TEXT_SIZE_((char *)" \n"));

    PrintStructElement(&SBGC32_Device, &RealTimeData.AxisRTD[ROLL].ACC_Data, "ACC Roll =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &RealTimeData.AxisRTD[PITCH].ACC_Data, "ACC Pitch =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &RealTimeData.AxisRTD[YAW].ACC_Data, "ACC Yaw =", _SIGNED_SHORT_);

    PrintStructElement(&SBGC32_Device, &RealTimeData.frameCamAngle[ROLL], "Roll Current Angle =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &RealTimeData.frameCamAngle[PITCH], "Pitch Current Angle =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &RealTimeData.frameCamAngle[YAW], "Yaw Current Angle =", _SIGNED_SHORT_);

    PrintStructElement(&SBGC32_Device, &RealTimeData.IMU_Temperature, "IMU Temperature =", _SIGNED_CHAR_);

    return SBGC32_Device._parserCurrentStatus;
}

TxRxStatus_t SBGC32_DemoControl(void)
{
    /* Getting adjvars values */
    /* Note: AdjVarGeneral.ID fields are already filled */
    SBGC32_GetAdjVarValues(&SBGC32_Device, AdjVarGeneral, countof_(AdjVarGeneral));

    /* Run the Demonstration Cycle */
    for (ui8 i = 0; i < 4; i++)
    {
        /* Printing. SBGC_ADJ_VARS_NAMES parameter must be SET_ON */
        for (ui8 k = 0; k < countof_(AdjVarGeneral); k++)
            PrintStructElement(&SBGC32_Device, &AdjVarGeneral[k].value, AdjVarGeneral[k].name, AdjVarGeneral[k].varType);

        Control.AxisC[YAW].angle = DEGREE_TO_ANGLE_INT(50);
        Control.AxisC[PITCH].angle = DEGREE_TO_ANGLE_INT(-20);
        SBGC32_Control(&SBGC32_Device, &Control);
        DELAY_MS_(5000);

        Control.AxisC[PITCH].angle = DEGREE_TO_ANGLE_INT(20);
        SBGC32_Control(&SBGC32_Device, &Control);
        DELAY_MS_(5000);

        Control.AxisC[YAW].angle = DEGREE_TO_ANGLE_INT(-50);
        SBGC32_Control(&SBGC32_Device, &Control);
        DELAY_MS_(5000);

        Control.AxisC[PITCH].angle = DEGREE_TO_ANGLE_INT(-20);
        SBGC32_Control(&SBGC32_Device, &Control);
        DELAY_MS_(5000);

        Control.AxisC[YAW].angle = DEGREE_TO_ANGLE_INT(0);
        Control.AxisC[PITCH].angle = DEGREE_TO_ANGLE_INT(0);
        SBGC32_Control(&SBGC32_Device, &Control);
        DELAY_MS_(5000);

        BeeperSettings.mode = BM_BEEPER_MODE_COMPLETE;
        SBGC32_PlayBeeper(&SBGC32_Device, &BeeperSettings);

        /* Adjustable Variables Re-Setting */
        for (ui8 k = 0; k < countof_(AdjVarGeneral); k++)
            /* Toggle Min : Max adjvars contrast */
            EditAdjVarValue(&AdjVarGeneral[k], ((i % 2 == 0) ? AdjVarGeneral[k].maxValue : AdjVarGeneral[k].minValue));

        SBGC32_SetAdjVarValues(&SBGC32_Device, AdjVarGeneral, countof_(AdjVarGeneral));
    }

    /* Saving all changed adjustable variables to EEPROM */
    /* SBGC32_SaveAllActiveAdjVarsToEEPROM(&SBGC32_Device);

    if (SBGC32_Device._confirmationParams.cmdID == CMD_SAVE_PARAMS_3)
        for (ui8 i = 0; i < countof_(AdjVarGeneral); i++)
            if (AdjVarGeneral[i].saveFlag != SAVED)
                AdjVarGeneral[i].saveFlag = SAVED; */

    /* or SBGC32_SaveAdjVarsToEEPROM(&SBGC32_Device, AdjVarGeneral, countof_(AdjVarGeneral)); */

    return SBGC32_Device._parserCurrentStatus;
}

void PrintDataStream(ui8 *pBuff)
{
    /* Preparing */
    ui8 BuffRPx = 2; // ui16 timestampMs offset

    BuffRPx += ConvertWithPM(RealTimeDataCustom.frameCamAngle, &pBuff[BuffRPx],
                             sizeof(RealTimeDataCustom.targetAngles), PM_DEFAULT_16BIT);
    BuffRPx += ConvertWithPM(RealTimeDataCustom.gyroData, &pBuff[BuffRPx],
                             sizeof(RealTimeDataCustom.gyroData), PM_DEFAULT_16BIT);
    BuffRPx += ConvertWithPM(RealTimeDataCustom.ACC_Data, &pBuff[BuffRPx],
                             sizeof(RealTimeDataCustom.ACC_Data), PM_DEFAULT_16BIT);

    /* Printing */
    PrintStructElement(&SBGC32_Device, &RealTimeDataCustom.frameCamAngle[ROLL], "Frame Camera Angle Roll =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &RealTimeDataCustom.frameCamAngle[PITCH], "Frame Camera Angle Pitch =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &RealTimeDataCustom.frameCamAngle[YAW], "Frame Camera Angle Yaw =", _SIGNED_SHORT_);

    PrintStructElement(&SBGC32_Device, &RealTimeDataCustom.gyroData[ROLL], "Gyro Roll =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &RealTimeDataCustom.gyroData[PITCH], "Gyro Pitch =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &RealTimeDataCustom.gyroData[YAW], "Gyro Yaw =", _SIGNED_SHORT_);

    PrintStructElement(&SBGC32_Device, &RealTimeDataCustom.ACC_Data[ROLL], "ACC Roll =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &RealTimeDataCustom.ACC_Data[PITCH], "ACC Pitch =", _SIGNED_SHORT_);
    PrintStructElement(&SBGC32_Device, &RealTimeDataCustom.ACC_Data[YAW], "ACC Yaw =", _SIGNED_SHORT_);

    PrintMessage(&SBGC32_Device, TEXT_SIZE_((char *)"__________________________\n\n"));
}