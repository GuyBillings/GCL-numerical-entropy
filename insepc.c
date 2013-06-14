/* C function calculate output separation
   Guy Billings UCL 2011*/
   
#include <mathlink.h>
#include <bpnfuncs.h>
///NEED TO MODIFY THE ARGUMENT LIST (SHORTEN)
extern double* insepc(int tc, int mu, int dini, int dinc, int dmax, double pini, double pinc, double pmax,int gamma, int prec);
double* insepc(int tc, int mu, int dini, int dinc, int dmax, double pini, double pinc, double pmax,int gamma,int prec)
{
   int i;
   int j;  
   int ddex;
   int pdex;
   //DIFFERENCE HERE IN HOW COMPILER INTERPRETS FOR LINUX AND MAC
   //int psiz=(int)((pmax-pini)/pinc)+1; //LINUX
   int psiz=round((pmax-pini)/pinc)+1;  //MAC
   
   mpz_t mump;
   mpz_init_set_ui(mump,(unsigned long int)mu);
 //  mpz_t dmp;
 //  mpz_init_set_ui(dmp,(unsigned long int)dmax);
 //  mpfr_t pinf;
 //  mpfr_init2(pinf,prec);//TEMP
 //  mpfr_set_d(pinf,0.4,MPFR_RNDN);
 //  mpz_t gamp;
 //  mpz_init(gamp);
 //  mpz_set_ui(gamp,gamma);
   //mpz_t ncodes;
   //mpz_init(ncodes);
   //mpz_set_ui(ncodes,tc);
  
   int pdfsize=(mu+1)*psiz;    
   mpfr_t * ipdfptr;
//   mpfr_t * ucodesptr;
//   mpfr_t * p0ptr;
   ipdfptr=(mpfr_t *) malloc(pdfsize*sizeof(mpfr_t));
   //ucodesptr=(mpfr_t *) malloc(pdfsize*sizeof(mpfr_t));
   //p0ptr=(mpfr_t *) malloc(pdfsize*sizeof(mpfr_t));
   for(i=0;i<=pdfsize-1;i++)
   {
    mpfr_init2(*(ipdfptr+i),prec);
    mpfr_set_d(*(ipdfptr+i),0,MPFR_RNDN);
   // mpfr_init2(*(ucodesptr+i),prec);
   // mpfr_init2(*(p0ptr+i),prec);
   }
   
 //  int rsize;
 //  double rtmp;
 //  if(dmax!=dini)   
 //  {
 //    rsize=(int)((mu+1)*(gamma+1)*psiz*((dmax-dini+1)*(dmax-dini+2)/2));
 //  }
 // else
 //  {
 //    rsize=(int)((mu+1)*(gamma+1)*psiz*dmax);
 //  }  
 //  mpfr_t * distptr;
 //  distptr=(mpfr_t *) malloc(rsize*sizeof(mpfr_t)); 
 //  for(i=0;i<=rsize-1;i++)
 //  {
 //    mpfr_init2(*(distptr+i),prec);
 //  }  
  
  // int esize;
  // if(dmax!=dini)
  // {
  //   esize=(int)psiz*((dmax-dini+1)*(dmax-dini+2)/2);
  // }
  // else
  // {
  //   esize=(int)psiz*dmax;
  // }
  
  // mpfr_t * ent;
  // ent=(mpfr_t *) malloc(esize*sizeof(mpfr_t)); 
  // for(i=0;i<=esize-1;i++)
  // {
  //   mpfr_init2(*(ent+i),prec);
  //   mpfr_set_d(*(ent+i),0,MPFR_RNDN);
  // }  
   
   //Binomial coeffs: (z,x) with z:0->gamma, x:0->gamma
  // mpfr_t * gamma_bcfr;
  // int gbf_size=(gamma+1)*(gamma+1);
  // gamma_bcfr=(mpfr_t *)malloc(gbf_size*sizeof(mpfr_t));
  // for(i=0;i<=gbf_size-1;i++)
  //   mpfr_init2(*(gamma_bcfr+i),prec);
        
   //Binomial coeffs: (y,x) with y:0->mu, x:0->d
  // int bacsize=(mu+1)*(dmax+1);
  // mpz_t beta_ai_cnt[bacsize];
  // for(i=0;i<=bacsize-1;i++)
  //   mpz_init(beta_ai_cnt[i]);   
     
   //Binomial coeffs: (y,x) with x:0->mu, y:0->mu
   mpfr_t mu_bcfr[(mu+1)*(mu+1)];
   for(i=0;i<=(mu+1)*(mu+1)-1;i++)
     mpfr_init2(mu_bcfr[i],prec);  
     
   //Binomial coeff: (mu,d) x:1->d (integer) 
   //mpz_t mu_bc[dmax]; 
   //for(i=0;i<=dmax-1;i++)  
   //  mpz_init(mu_bc[i]);
     
   //int ppsize=0.5*(dmax*(dmax+1)-(dini-1)*dini)*(mu+1);
   //mpfr_t * ppptr;
   //ppptr=(mpfr_t *) malloc(ppsize*sizeof(mpfr_t)); 
   //for(i=0;i<=ppsize-1;i++)
   //{
   //  mpfr_init2(*(ppptr+i),prec);  
   //  mpfr_set_d(*(ppptr+i),0,MPFR_RNDN);
   //} 
   
   int psepsize=psiz;
   
   mpfr_t *psep;
   psep=(mpfr_t *) malloc(psepsize*sizeof(mpfr_t));
   for(i=0;i<=psepsize-1;i++)
   {
     mpfr_init2(*(psep+i),prec);  
     mpfr_set_d(*(psep+i),0,MPFR_RNDN);
   } 
     
   mlextended_double *retvals;
   retvals=(mlextended_double *) malloc(psepsize*sizeof(mlextended_double)); 
   
//   tabulate_bins_fr(gamma_bcfr,0,gamma,0,gamma,prec);
   tabulate_bins_fr(mu_bcfr,0,mu,0,mu,prec);
//   tabulate_bins_z(mu_bc,mu,mu,1,dmax);           
//   tabulate_bins_z(beta_ai_cnt,0,mu,0,dmax);         
   inpdf(ipdfptr,mump,psiz,pini,pinc,mu,(mu_bcfr+mu*(mu+1)),prec);
//   uniquecodes(ucodesptr,psiz,mump,mu,ncodes,ipdfptr,prec);
//   pphi(ppptr,mump,dmp,mu,dini,dmax,beta_ai_cnt,mu_bc,prec);
//   scjoint(distptr,mu,gamma,psiz,dini,dmax,(gamma_bcfr+gamma*(gamma+1)),ipdfptr,ppptr,prec);
      
//   for(i=0;i<=pdfsize-1;i++)
//   {
//      if(mpfr_cmp_si(*(ucodesptr+i),0)<=0)
//        mpfr_set_ui(*(p0ptr+i),0,MPFR_RNDN);       
//      else 
//        mpfr_ui_div(*(p0ptr+i),1,*(ucodesptr+i),MPFR_RNDN);
//   }  
   
   patt_sep_in(psep,ipdfptr,mu_bcfr,psiz,mu,pinc,prec);
   
   for(i=0;i<=psepsize-1;i++)
   {
     *(retvals+i)=(mlextended_double)mpfr_get_d(*(psep+i),MPFR_RNDN);
   } 
   MLPutReal128List(stdlink, (mlextended_double *)retvals,psepsize);
   
   mpz_clear(mump);
   //mpz_clear(dmp);
   //mpfr_clear(pinf);
   //mpz_clear(gamp);
   //mpz_clear(ncodes);
   
   for(i=0;i<=pdfsize-1;i++)
   {
    mpfr_clear(*(ipdfptr+i));
 //   mpfr_clear(*(ucodesptr+i));
 //   mpfr_clear(*(p0ptr+i));
   }
   
   //for(i=0;i<=rsize-1;i++)
   //  mpfr_clear(*(distptr+i));

   for(i=0;i<=psepsize-1;i++)
     mpfr_clear(*(psep+i));

  //for(i=0;i<=(gamma+1)*(gamma+1)-1;i++)
  //   mpfr_clear(*(gamma_bcfr+i));
     
  // for(i=0;i<=bacsize-1;i++)
  //   mpz_init(beta_ai_cnt[i]);     
     
   for(i=0;i<=(mu+1)*(mu+1)-1;i++)
     mpfr_clear(mu_bcfr[i]);    
        
  // for(i=0;i<=dmax-1;i++)  
  //   mpz_clear(mu_bc[i]);
     
  // for(i=0;i<=ppsize-1;i++)
  //   mpfr_clear(*(ppptr+i));  
  
   return 0;
}

int main(int argc, char* argv[])
{
	return MLMain(argc, argv);
}

