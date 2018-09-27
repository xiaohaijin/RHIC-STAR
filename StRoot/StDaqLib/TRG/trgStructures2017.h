#ifndef trgStructures2017_h
#define trgStructures2017_h
/******
*
*     Layout of new Trigger Data Block
*
*     J.M. Nelson        30 January 2009
*
*     Notes:  The event descriptor will describe data from either
*     the Mk1 or Mk2 TCUs.   The variable TCU_Mark will be 1 for Mk1
*     and 2 for the Mk2 TCU.  Variables not used by one or other of the
*     TCUs will be zero.   
*
*     The data block structure will always begin with a 4 character
*     name, followed by the byte-count of data following.  The structure of
*     data will depend on the configuration of particular crates.  
*
*     Note:  PrePost data will only be available on local trigger disks and
*     will not be present in event files.
******************************************************************************/
#define y17FORMAT_VERSION        0x16120844      /* Format: yymmddvv */
#define y17MAX_TRG_BLK_SIZE          122896      /* Current total: 113.25k bytes for pre/post non-zero suppressed data.  Allow 120k */
#define y17MAX_OFFLEN                    20      /* Depends on the number of crates in the system */

#define   y17L1_CONF_NUM       1
#define  y17BC1_CONF_NUM       2
#define  y17MXQ_CONF_NUM       3
#define  y17MIX_CONF_NUM       4
#define  y17BCW_CONF_NUM       5
#define  y17BCE_CONF_NUM       6
#define  y17EPQ_CONF_NUM       7
#define  y17BBC_CONF_NUM       8
#define  y17BBQ_CONF_NUM       9
#define  y17FMS_CONF_NUM      10
#define  y17QT1_CONF_NUM      11
#define  y17QT2_CONF_NUM      12
#define  y17QT3_CONF_NUM      13
#define  y17QT4_CONF_NUM      14
#define  y17FQ1_CONF_NUM      15
#define  y17FQ2_CONF_NUM      16

#define y17ADD_BIT_FORCE          5              /* Force store of this event */
#define y17ADD_BIT_L2_5           6              /* Level 2.5 abort */
#define y17ADD_BIT_SIM            7              /* Simulated event - used by DAQ */

#define L2RESULTS_2017_OFFSET_EMC_CHECK   1
#define L2RESULTS_2017_OFFSET_EMC_PED     2
#define L2RESULTS_2017_OFFSET_BGAMMA      3
#define L2RESULTS_2017_OFFSET_EGAMMA      6
#define L2RESULTS_2017_OFFSET_DIJET       9 
#define L2RESULTS_2017_OFFSET_UPSILON     17
#define L2RESULTS_2017_OFFSET_BEMCW       20
#define L2RESULTS_2017_OFFSET_BHIEN       42
#define L2RESULTS_2017_OFFSET_EHIEN       0
#define L2RESULTS_2017_OFFSET_BTOW_CAL    0
#define L2RESULTS_2017_OFFSET_ETOW_CAL    0

    /* Event Descriptor Data Structures */
    
//#pragma pack(1)

typedef struct {
  char           name[3];                     /* Contains  EVD */
  char           TrgDataFmtVer;               /* Exception for use by DAQ (LS byte of FORMAT_VERSION) */
  int            length;                      /* Byte count of data that follows */
  unsigned int   bunchXing_hi;
  unsigned int   bunchXing_lo;                /* Two parts of RHIC bunch crossing number */
  unsigned short actionWdDetectorBitMask;     /* from Fifo 1 */
  unsigned char  actionWdTrgCommand;          /* from Fifo 1 */
  unsigned char  actionWdDaqCommand;          /* from Fifo 1 */  
  unsigned short TrgToken;                    /* from Fifo 2 */
  unsigned short addBits;                     /* used by trigger/daq: bit 5=Force store; bit 6=L2.5 abort; bit 7=1 is fake data */
  unsigned short DSMInput;                    /* only for use with Mk1 TCU.  0 if Mk2 TCU is used */
  unsigned short externalBusy;                /* from Fifo 9 (Fifo 3 Mk1 TCU) */
  unsigned short internalBusy;                /* from Fifo 9 (Mk2 TCU) */

  unsigned short trgDetMask;                  /* was physicsWord */
  unsigned short tcuCtrBunch_hi;              /* was TriggerWord */

  unsigned short DSMAddress;                  /* from Fifo 10 (Fifo 6 Mk1 TCU) */
  unsigned short TCU_Mark;                    /* TCU_Mark Mk1=1 Mk2=2 */
  unsigned short npre;                        /* (crate_mask & 0xfff) << 4 | npre  */
  unsigned short npost;                       /* (crate_mask & 0xfff000)>>8| npost */
  unsigned short res1;                        /* (crate_mask & 0xff000000)>>20 | res1&0xf */
} EvtDescData2017;

//#pragma pack()

      /* L1 DSM data structures */

typedef struct {
  char               name[4];                 /* Contains  L1DS */
  int                length;                  /* Byte count of data that follows */
  unsigned short     TOF[8];                  /* TOF and MTD data */
  unsigned short     VTX[8];                  /* Separate VPD, ZDC and BBC DSMs have been replaced with this one */
  unsigned short     EMC[8];                  /* Contents of 1 EMC IB - results of separate BEMC and EEMC DSMs */
  unsigned short     TPCMask[8];              /* TPC mask for DAQ10K */
  unsigned short     BCdata[16];              /* Contents of 2 Bunch Crossing DSMs IB's */       
  unsigned short     specialTriggers[8];      /* Contents of 1 Special Trigger DSM - all the special trigger requests */
  unsigned short     FPD[8];                  /* Contents of 1 FMS and FPD IB */
  unsigned short     lastDSM[8];              /* Contents of last DSM IB - results of all DSM trees */
} L1_DSM_Data2017;

      /* Trigger Summary Data Structures */

typedef struct {
  char           name[4];                     /* Contains  TSUM */
  int            length;                      /* Byte count of data that follows */
  unsigned int   L1Sum[2];                    /* L1 Summary */
  unsigned int   L2Sum[2];                    /* L2 Summary */
  unsigned int   L1Result[32];                /* Result from L1 CPU */
  unsigned int   L2Result[64];                /* Result from L2 CPU */
  unsigned int   C2Result[64];                /* Result from last algorithm */
  unsigned int   LocalClocks[32];	      /* localClock values from RCC2*/
} TrgSumData2017;

typedef struct {
  char name[4];
  int length;                                 /* Byte count of data that follows */
  unsigned int data[1];                       /* NB: this definition is generic but would vary depending on actual data */
} DataBlock2017;

typedef struct {
  char name[4];                               /* Contains BBC */
  int length;                                 /* Byte count of data that follows */
  unsigned short BBClayer1[16];               /* This is the layer1 DSM that feeds the VTX DSM + EPD */
  unsigned short ZDClayer1[8];                /* This is the new layer1 ZDC DSM that also feeds the VTX DSM */
  unsigned short VPD[8];                      /* ADC & TAC values for VPD detectors*/
} BBCBlock2017;

typedef struct {
  char name[4];                               /* Contains MIX */
  int length;                                 /* Byte count of data that follows */
  unsigned short FPDEastNSLayer1[8];          /* FPD east north/south layer 1  -> pp2pp */  
  unsigned char  MTD_P2PLayer1[16];           /* Data from MTD and PP2PP */
  unsigned short TOFLayer1[8];                /* This is TOF Layer 1 */
  unsigned short TOF[48];                     /* TOF data */
  unsigned short TPCpreMask[24];              /* EMC, MTD, & TOF TPC Grid Masks -> HCAL */
} MIXBlock2017;

typedef struct  {
  char name[4];
  int length;                                 /* Byte count of data that follows */
  int dataLoss;                               /* Byte count of data truncated due to buffer limitations */
  unsigned int data[1];                       /* NB: this definition is generic but would vary depending on actual data */
} QTBlock2017;

typedef struct {
  char name[4];
  int length;
  unsigned char BEMCEast[240];                /* 15 DSMs covering the East half of BEMC */
} BEastBlock2017;

typedef struct {
  char name[4];
  int length;
  unsigned char BEMCWest[240];                /* 15 DSMs covering the West half of BEMC */
} BWestBlock2017;

typedef struct {
  char name[4];
  int length;
  unsigned short BEMClayer1[48];              /* 6 DSMs for BEMC at layer1 */
  unsigned short EEMClayer1[16];              /* 2 DSMs for EEMC at layer1 */
  unsigned char  EEMC[144];                   /* 9 DSMs for EEMC at layer0 */
} BELayerBlock2017;

typedef struct {
  char name[4];
  int length;
  unsigned char FMS[256];                     /* 16 DSMs for FMS */
} FMSBlock2017;

typedef struct {
  int offset;                                 /* Offset (in bytes) from the start of Trigger block to data */
  int length;                                 /* Length (in bytes) of data */
} TrgOfflen2017; 

typedef struct {
  int FormatVersion;                          /* Trigger Data Definition Version yymmddvv */
  int totalTriggerLength;                     /* Total length (bytes) of complete Trigger Block */
  int eventNumber;                            /* Event number in this run */
  TrgOfflen2017 EventDesc_ofl;                    /* Offset/length pair to Event Descriptor */
  TrgOfflen2017 L1_DSM_ofl;                       /* Offset/length pair to L1 DSM Data */
  TrgOfflen2017 Summary_ofl;                      /* Offset/length pair to Summary Data */
  TrgOfflen2017 MainX[y17MAX_OFFLEN];                /* Offset/length pairs for main crossing */
  int PrePostList[10];                        /* Offsets to offset/length pairs to Pre and Post crossing */
  int raw_data[y17MAX_TRG_BLK_SIZE/4];           /* Storage for raw data */
} TriggerDataBlk2017;

#endif

