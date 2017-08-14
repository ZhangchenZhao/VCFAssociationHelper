#ifndef INT8_MAX
#define 	INT8_MAX   0x7f
#endif 

#ifndef INT8_MIN
#define 	INT8_MIN   (-INT8_MAX - 1)
#endif 

#ifndef UINT8_MAX   
#define 	UINT8_MAX   (INT8_MAX * 2 + 1)
#endif 

#ifndef  INT16_MAX   
#define 	INT16_MAX   0x7fff
#endif 

#ifndef  INT16_MIN   
#define 	INT16_MIN   (-INT16_MAX - 1)
#endif 
 
#ifndef UINT16_MAX   
#define 	UINT16_MAX   (__CONCAT(INT16_MAX, U) * 2U + 1U)
#endif 
 
#ifndef INT32_MAX   
#define 	INT32_MAX   0x7fffffffL
#endif 
 
#ifndef INT32_MIN   
#define 	INT32_MIN   (-INT32_MAX - 1L)
#endif 
 
#ifndef UINT32_MAX   
#define 	UINT32_MAX   (__CONCAT(INT32_MAX, U) * 2UL + 1UL)
#endif 
 
#ifndef INT64_MAX   
#define 	INT64_MAX   0x7fffffffffffffffLL
#endif 
 
#ifndef INT64_MIN   
#define 	INT64_MIN   (-INT64_MAX - 1LL)
#endif 
 
#ifndef UINT64_MAX   
#define 	UINT64_MAX   (__CONCAT(INT64_MAX, U) * 2ULL + 1ULL)
#endif 
