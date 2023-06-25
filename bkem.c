#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "bkem.h"
#include <time.h>
clock_t setup_time, keygen_time, enc_time, dec_time;
void setup_global_system(bkem_global_params_t *gps, const char *pstr, int N) {
    
    bkem_global_params_t params;
    params = pbc_malloc(sizeof(struct bkem_global_params_s));

    params->N = N;
    
    pairing_init_set_str(params->pairing, pstr);

    *gps = params;
}

/*void generateRandomString(int length, int* randomArray) {
    srand(time(NULL));
    
    for (int i = 0; i < length; i++) {
        randomArray[i] = rand() % 2;
    }
}
*/
void generateRandomArrays(int numArrays, int arrayLength, int randomArrays[Max_N][LogMax_N]) {
    srand(time(NULL));
    
    for (int i = 0; i < numArrays; i++) {
        for (int j = 0; j < arrayLength; j++) {
            randomArrays[i][j] = rand() % 2;
        }
    }
}


int sizeOf(ID ij)
{
    return (int)(2*sizeof(int));
}

void hashID(element_t hash, ID IDij, bkem_global_params_t gps)
{	
    element_init_Zr(hash, gps->pairing);
    element_from_hash(hash, &IDij, sizeOf(IDij));
}


void setup(bkem_system_t *sys, bkem_global_params_t gps) 
{
    setup_time = clock();
    
    // 
    bkem_system_t gbs;
    bkem_secret_key_t sk;
    gbs = pbc_malloc(sizeof(struct bkem_system_s));
    gbs->PK = pbc_malloc(sizeof(struct pubkey_s));

    // ---------------------------------Choose random generator g --------------------------------------------
    element_init_G1(gbs->PK->g, gps->pairing);
    element_random(gbs->PK->g);

    //----------------------------------Choose another generator ghat=gg--------------------------------------
    element_init_G2(gbs->PK->gg, gps->pairing);
    element_random(gbs->PK->gg);

    // ---------------------------------random alpha in Zn ---------------------------------------------------
    element_t alpha;
    element_init_Zr(alpha, gps->pairing);
    element_random(alpha);

    // ---------------------------------random beta in Zn-----------------------------------------------------
    element_t beta;
    element_init_Zr(beta, gps->pairing);
    element_random(beta);

    // ---------------------------------random xhat in Zn ----------------------------------------------------
    element_t xhat;
    element_init_Zr(xhat, gps->pairing);
    element_random(xhat);

    // ---------------------------------random yhat in Zn ----------------------------------------------------
    element_t yhat;
    element_init_Zr(yhat, gps->pairing);
    element_random(yhat);
    
    

    /*
    element_printf("alpha = %B\n", alpha);
    element_printf("beta = %B\n", beta);
    element_printf("xhat = %B\n", xhat);
    element_printf("yhat = %B\n\n", yhat);
    */	

   // -------------------------------Compute the component of MPK ---------------------------------------------
   gbs->PK->mpk_i = pbc_malloc( 5 * sizeof(element_t));
   int size_of_MPK=(Max_N+LogMax_N+5) * sizeof(element_t);		// 
   //element_printf("size_of_MPK = %d in bytes\n\n", size_of_MPK);
    
   // element_printf("Compute the component of MPK\n");
	
    //-----------------------------Set the first element to g--------------------------------------------------
    element_init_G1(gbs->PK->mpk_i[0], gps->pairing);
    element_set(gbs->PK->mpk_i[0],gbs->PK->g);
    //element_printf("g = %B\n\n", gbs->PK->mpk_i[0]);
    
    
    //-------------------------------Set the first element to Omega=e(g,ghat)^{alpha.beta}--------------------------------
    element_t t1,t2,t3;   
    element_init_GT(gbs->PK->mpk_i[1], gps->pairing);
    element_init_Zr(t1, gps->pairing);
    element_mul(t1,alpha,beta);
    pairing_apply(gbs->PK->mpk_i[1], gbs->PK->g, gbs->PK->gg,gps->pairing);
    element_pow_zn(gbs->PK->mpk_i[1], gbs->PK->mpk_i[1], t1);
    //element_printf("Omega = %B\n\n", gbs->PK->mpk_i[1]);
    
    //-------------------------------Set the first element to Lambda=e(g,ghat)^{alpha.(beta-1)}--------------------------------  
    element_init_GT(gbs->PK->mpk_i[2], gps->pairing);
    element_init_Zr(t2, gps->pairing);
    element_init_Zr(t3, gps->pairing);
    element_mul(t2,alpha,beta);
    element_sub(t3,t2,alpha);
    pairing_apply(gbs->PK->mpk_i[2], gbs->PK->g, gbs->PK->gg,gps->pairing);
    element_pow_zn(gbs->PK->mpk_i[2], gbs->PK->mpk_i[2], t3);
    //element_printf("Lambda = %B\n\n", gbs->PK->mpk_i[2]);
    
    //-------------------------------Set the first element to Xhat=(ghat)^xhat-------------------------------------
    element_init_G2(gbs->PK->mpk_i[3], gps->pairing);
    element_pow_zn(gbs->PK->mpk_i[3], gbs->PK->gg, xhat);
    //element_printf("Xhat = %B\n\n", gbs->PK->mpk_i[3]);
    
    //------------------------------Set the first element to Yhat=(ghat)^yhat----------------------------------------
    element_init_G2(gbs->PK->mpk_i[4], gps->pairing);
    element_pow_zn(gbs->PK->mpk_i[4], gbs->PK->gg, yhat);
    //element_printf("Yhat = %B\n\n", gbs->PK->mpk_i[4]);
    
    //---------------------------------Computes T_i and TT_i------------------------------------------------------------------
    for (int j = 0; j < LogMax_N; ++j) 
    {
	element_init_Zr(gbs->PK->t[j], gps->pairing);
    	element_random(gbs->PK->t[j]);
    	//element_printf("t[%d] = %B\n\n", j, gbs->PK->t[j]);
    	element_init_G1(gbs->PK->T[j], gps->pairing);
    	element_init_G2(gbs->PK->TT[j], gps->pairing);
    	element_pow_zn(gbs->PK->T[j], gbs->PK->g, gbs->PK->t[j]);
    	element_pow_zn(gbs->PK->TT[j], gbs->PK->gg, gbs->PK->t[j]);
  	//element_printf("T[%d] = %B\n\n", j, gbs->PK->T[j]);
     	//element_printf("TT[%d] = %B\n\n", j, gbs->PK->T[j]);
    }
    
    
   //================================Master Secret-key==============================================================
    //element_printf("Compute the component of MSK\n");
    gbs->PK->msk_i = pbc_malloc(2 * sizeof(element_t));
    
    //-------------------------------Set the first element to ghat---------------------------------------------
    element_init_G2(gbs->PK->msk_i[0], gps->pairing);
    element_set(gbs->PK->msk_i[0],gbs->PK->gg);
    //element_printf("ghat = %B\n\n", gbs->PK->msk_i[0]);
    
    //-------------------------------Set the first element to ghat_1=(ghat)^alpha--------------------------------------
    element_init_G2(gbs->PK->msk_i[1], gps->pairing);
    element_pow_zn(gbs->PK->msk_i[1], gbs->PK->gg, alpha);
    //element_printf("ghat_1 = %B\n\n", gbs->PK->msk_i[1]);
   
/*    setup_time = clock() - setup_time;
    double time_taken0 = ((double)setup_time)/CLOCKS_PER_SEC; // in seconds 
    printf("Setup took %f seconds to execute \n\n", time_taken0);  
 */   
     
    //int size_of_MSK=3 * sizeof(element_t);
    //element_printf("size_of_MSK = %d in bytes\n\n", size_of_MSK);
   
   //------------------MPK and MSK generation is done ----------------------------------------------------------------
  
 	

// -----------To Compute private keys of users--------------------------------------------------------

    keygen_time = clock();
    int user=25;
    element_t d, d1,e1;
/*    int* randomArray = (int*)malloc(LogMax_N * sizeof(int)); // Alocation for hash value
    generateRandomString(LogMax_N, randomArray);
    printf("The hash value of %d-user: ", user);
    for (int j = 0; j < LogMax_N; j++) {
        printf("%d", randomArray[j]);
    }
    printf("\n");
*/   
    	 

    int hash_value[Max_N][LogMax_N];
    generateRandomArrays(Max_N, LogMax_N, hash_value);
/*    for (int i = 0; i < Max_N; i++) {
    	printf("The hash value of %d-user: ", i);
    	for (int j = 0; j < LogMax_N; j++) {
        	printf("%d", hash_value[i][j]);
    	}
    	printf("\n");
    }
*/   
    
 
    for(int i=0;i<Max_N;i++){
    element_init_G2(gbs->Hat_FF[i], gps->pairing);
    element_set1(gbs->Hat_FF[i]); 		
    for(int j=0;j<LogMax_N;j++)
    {
    	if(hash_value[i][j]==1)
    	{
    		//element_printf("z[%d] = %B\n\n", j,gbs->PK->z[j]);
    		element_mul(gbs->Hat_FF[i],gbs->Hat_FF[i],gbs->PK->TT[j]);
    	}    		
    }
    //element_printf("gbs->Hat_FF[%d] = %B\n\n", i, gbs->Hat_FF[i]);
       
   }
    
  //------------------------ Compute the private keys SK_j_i -----------------------------------------------------
    
    for(int i=0;i<Max_N;i++){
    element_t r_user, Hat_F_user1,sk_1,sk_0; 
    element_init_G2(gbs->SK[i][0], gps->pairing);
    element_init_Zr(gbs->r[i], gps->pairing);
    element_random(gbs->r[i]);
    element_init_G2(Hat_F_user1, gps->pairing);
    element_pow_zn(Hat_F_user1,gbs->Hat_FF[i], gbs->r[i]);
    element_mul(gbs->SK[i][0],gbs->PK->msk_i[1], Hat_F_user1);
    //element_printf("The 1st component of secret-key d[%d][0] = %B\n\n", i, gbs->SK[i][0]);

   
    
    element_init_G2(gbs->SK[i][1], gps->pairing);
    element_pow_zn(gbs->SK[i][1],gbs->PK->msk_i[0], gbs->r[i]);
    //element_printf("The 2nd component of secret-key d[%d][1] = %B\n\n", i, gbs->SK[i][1]);
    }
   
 /*   keygen_time = clock() - keygen_time;
    double time_taken1 = ((double)keygen_time)/CLOCKS_PER_SEC; // in seconds 
    printf("KeyGen took %f seconds to execute \n\n", time_taken1); 
 */

    *sys = gbs;
    element_clear(alpha);
    element_clear(beta);
    element_clear(xhat);
    element_clear(yhat);  
    
    //----------------------------Key Gen is done ----------------
}



//=================================Encryption =============================================




void get_enc_key( bkem_system_t gbs, bkem_global_params_t gps) 
{
     enc_time=clock();
     element_t s,temp1,temp2,temp3,temp4;
     element_init_Zr(s, gps->pairing);
     element_random(s);//element_printf("C_0= %B\n\n", gbs->C_0);
     element_init_GT(gbs->M, gps->pairing);
     element_random(gbs->M);		//Assign the broadcasting message 
     //element_printf("The message which I want to broadcast= %B\n\n", gbs->M);
     element_init_GT(temp1, gps->pairing);
     element_pow_zn(temp1, gbs->PK->mpk_i[1], s);
     element_init_GT(gbs->C_0, gps->pairing);		// Initialize the 1st ciphertext component
     element_mul(gbs->C_0, gbs->M, temp1);		//Calculate the 1st ciphertext component
     //element_printf("C_0= %B\n\n", gbs->C_0);
     element_init_GT(gbs->C_1, gps->pairing);		// Initialize the 2nd ciphertext component
     element_pow_zn(gbs->C_1, gbs->PK->mpk_i[2], s);		//Calculate the 2nd ciphertext component
     //element_printf("C_1= %B\n\n", gbs->C_1);
     element_init_G1(gbs->C_2, gps->pairing);		// Initialize the 2nd ciphertext component
     element_pow_zn(gbs->C_2, gbs->PK->g, s);		//Calculate the 2nd ciphertext component
     //element_printf("C_2= %B\n\n", gbs->C_2);
     
    int hash_value[Max_N][LogMax_N];
    generateRandomArrays(Max_N, LogMax_N, hash_value);
/*    for (int i = 0; i < Subs_Num; i++) {
    	printf("The hash value of %d-user: ", i);
    	for (int j = 0; j < LogMax_N; j++) {
        	printf("%d", hash_value[i][j]);
    	}
    	printf("\n");
    }
*/  
   		
    for(int i=0;i<Subs_Num; i++){
    	element_init_G1(gbs->Hat_F[i], gps->pairing);
    	element_set1(gbs->Hat_F[i]);
    	for(int j=0;j<LogMax_N;j++)
    	{
    		if(hash_value[i][j]==1)
    		{
    		//element_printf("z[%d] = %B\n\n", j,gbs->PK->z[j]);
    		element_mul(gbs->Hat_F[i],gbs->Hat_F[i],gbs->PK->T[j]);
    		}    		
    	}	
    	//element_printf("Hat_F[%d] = %B\n\n", i,gbs->Hat_F[i]);
    	element_init_G1(gbs->C_3[i], gps->pairing);
    	element_pow_zn(gbs->C_3[i],gbs->Hat_F[i],s);
    	//element_printf("The 3rd component C_3[%d] = %B\n\n", i,gbs->C_3[i]);
    }
    
    	element_init_Zr(gbs->theta, gps->pairing);
    	element_random(gbs->theta);
    	element_init_G2(temp2, gps->pairing);
    	element_pow_zn(temp2,gbs->PK->mpk_i[4],gbs->theta);
    	element_mul(temp2,temp2,gbs->PK->mpk_i[3]);
    	element_init_G2(gbs->hat_Gamma, gps->pairing);
    	element_pow_zn(gbs->hat_Gamma,temp2,s);
    	
    	
    	
    
    	enc_time = clock() - enc_time; 
	double time_taken2 = ((double)enc_time)/(CLOCKS_PER_SEC); // in seconds 
  	printf("Encryption algorithm took %f seconds to execute for %d users \n\n", time_taken2, Subs_Num); 
    
       element_clear(s);         	
}


void get_decryption_key(bkem_global_params_t gps, bkem_system_t gbs, pubkey_t PK)
 {
 	dec_time=clock();
 	int user=7;
 	//element_printf("The 1st component of secret-key d[%d][0] = %B\n\n", user, gbs->SK[user][0]);
 	//element_printf("The 1st component of secret-key d[%d][0] = %B\n\n", user, gbs->SK[user][1]);
 		
	element_t t1,t2,t3,t4,t5,t6,Guess_M;
 	element_init_GT(t1, gps->pairing);
 	element_div(t1,gbs->C_0, gbs->C_1);
 	element_init_GT(t2, gps->pairing); 	
 	pairing_apply(t2,gbs->C_3[user], gbs->SK[user][1],gps->pairing);
 	element_init_GT(t3, gps->pairing);
     	pairing_apply(t3,gbs->C_2, gbs->SK[user][0],gps->pairing);
     	element_init_GT(t4, gps->pairing);
     	element_div(t4, t2, t3);
     	element_init_GT(Guess_M, gps->pairing); 	
 	element_mul(Guess_M,t1, t4); 
 		
 	
 	dec_time = clock() - dec_time; 
   	double time_taken3 = ((double)dec_time)/(CLOCKS_PER_SEC); // in seconds 
  	printf("Decryption algorithm took %f seconds to execute for a subscribed user\n\n", time_taken3);
  	
  	//element_printf("The plaintext is %B\n\n",gbs->M);
  	//element_printf("The recover message is %B\n\n",gbs->M); 
 	//element_printf("The recover message is %B\n\n",Guess_M); 
 
}


/*void free_global_params(bkem_global_params_t gbs) {
    if (!gbs)
        return;

    pairing_clear(gbs->pairing);
    free(gbs);
}
*/
/*void free_pubkey(pubkey_t pk, bkem_global_params_t gbs) {
    if (!pk)
        return;

    element_clear(pk->g);

    int i;
    for (i = 0; i <= gbs->N; ++i) {
        element_clear(pk->g_i[i]);
    }

    //for (i = 0; i < gbs->A; ++i) {
       // element_clear(pk->v_i[0]);
    //}

}
*/
/*void free_bkem_system(bkem_system_t sys, bkem_global_params_t gbs) {
    if (!sys)
        return;

    free_pubkey(sys->PK, gbs);

    int i;
    for (i = 0; i < gbs->N; ++i) {
        element_clear(sys->d_i[i]);
    }*/
//}
