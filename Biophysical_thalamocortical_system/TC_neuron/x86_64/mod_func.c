#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _HH2_reg(void);
extern void _IT_reg(void);
extern void _cadecay_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," HH2.mod");
    fprintf(stderr," IT.mod");
    fprintf(stderr," cadecay.mod");
    fprintf(stderr, "\n");
  }
  _HH2_reg();
  _IT_reg();
  _cadecay_reg();
}
