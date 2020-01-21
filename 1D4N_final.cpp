#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <mpi.h>
////////////////////
#define n 1000               ///Ввод размерности решетки
float Polevar = 0;
//float Polevar = 0.01;           // Поле
/*long int prohod_MC = 1e8; /// Ввод числа шагов Монте-Карло
unsigned long long between_exchange=1e6; //между репличными обменами
*/
unsigned long long  Prohod_MC_equilibration = 1e8;      /// Ввод числа шагов разогрева
unsigned long long  Prohod_MC_sampling = 1e8;       /// Ввод числа шагов Монте-Карло
unsigned long long between_exchange_equilibration=1e3; //между репличными обменами
unsigned long long between_exchange=1e3; //между репличными обменами
bool replica_on=1;
int current_value = 0;    // выводить текущие значения или нет
unsigned long long CountOut = 1e8;  // количество вывода

///////////////////////

using namespace std;

///////////////////////////////////     Ввод переменных
int num_podhod=10;                      // количвество попыток найти систему связей с мин. кол-м ошибок.
double mintemp=0.1;                     // выбор min температуры
double maxtemp=3;                       // выбор max температуры
float Dt = 0.1;                         // шаг температуры
//int CountOut = 1e6;                // количество вывода



////////////////////////////////////
void Slspin (short *spin);
void SpinFerr(short *spin);
void Vivodspin(int spin);
int Slsvaz( short *svaz, int num_podhod);
void VivodSSS( short *spin,short *svaz);
void Energy(short *spin,short *svaz,float *E);
float EnergyIJ(int i, short *spin, short *svaz);
void VivodEnergy( float *E);
void CopyE(float *E,float *E1);
int Claster( float *E1);
int ClasterFerr( float *E1);
void ColorOut2( float *E, short *spin);
void ColorOut( float *E, short *spin);
float SumEnergy(float *E);
void ColorOutIJ( float *E, short *spin, int i,int j);
void ColorOutSpin( float *E, short *spin);
int MaxClass(int per, float *E1, int tmp, int *Ochered);
int MaxClassFerr(int per, float *E1, int tmp, int *Ochered);
int sosed_sl(int top);
int sosed_sp(int top);
int sosed_spsp(int top);
int sosed_slsl(int top);
double SROTKL(double OutputMass, double SRKVOTKL);
double MagnMet( short *spin);
void PhazDiagramm(int var04, int var02, int var00, int var20, int var40, float *E);
void MCL(double mintemp, double maxtemp, short *spin, short *svaz, float *E, float *E1, int NUM_step_mcl);
void MCProhod(double temp,short *spin,short *svaz,float *E);
void MCProhod2(double temp,short *spin,short *svaz,float *E, int rank);

void Slspin(short *spin)
{
    int v=0;
    for(int i=0;i<n;i++)
    {
        v=rand()%(2);
        if(v)
            spin[i]=1;
        else
            spin[i]=-1;
    }
}
//////////////////////////////////////////////////// Забивает матрицу спинов 1
void SpinFerr(short *spin)
{
    for(int i=0;i<n;i++)
    {
        spin[i]=1;
    }
}

///////////////////////////////////////////////////////Подсчитывает энегрию системы
void Energy(short *spin,short *svaz,float *E)
{
    for (int i = 2; i < (n-2); i++)
    {
        E[i] = -svaz[i] * spin[i] * spin[i+1] - svaz[i-1] * spin[i] * spin[i-1]- svaz[n+i] * spin[i] * spin[i+2] - svaz[(n+i)-2] * spin[i] * spin[i-2] - Polevar;

    }
        E[0] = -svaz[0] * spin[0] * spin[1] - svaz[n-1] * spin[0] * spin[n-1] - svaz[n] * spin[0] * spin[2] - svaz[2*n-2] * spin[0] * spin[n-2] - Polevar;
        E[1] = -svaz[0] * spin[0] * spin[1] - svaz[1] * spin[1] * spin[2] - svaz[n+1] * spin[1] * spin[3] - svaz[2*n-1] * spin[1] * spin[n-1] -  Polevar;
        E[n-1] = -svaz[n-1] * spin[n-1] * spin[0] - svaz[n-2] * spin[n-1] * spin[n-2] - svaz[2*n-1] * spin[n-1] * spin[1] - svaz[n-1] * spin[n-1] * spin[0] - Polevar;
        E[n-2] = -svaz[n-2] * spin[n-2] * spin[n-1] - svaz[n-3] * spin[n-3] * spin[n-2] - svaz[2*n-2] * spin[n-2] * spin[0] - svaz[2*n-4] * spin[n-2] * spin[n-4] - Polevar;


}
///////////////////////////////////////////////////////Копирует массив с энергией Е в Е1
void CopyE(float *E,float *E1)
{
    for (int i = 0; i < n; i++)
    {
        E1[i] = E[i];

    }
}
///////////////////////////////////////////////////////Копирует массив с энергией Е в Е1
void CopyE(short *E,float *E1)
{
    for (int i = 0; i < n; i++)
    {
        E1[i] = E[i];

    }
}
///////////////////////////////////////////////////////Копирует массив с энергией Е в Е1
void CopyE(float *E,short *E1)
{
    for (int i = 0; i < n; i++)
    {
        E1[i] = E[i];

    }
}
///////////////////////////////////////////////////////Копирует массив с энергией Е в Е1
void CopyE(short *E,short *E1)
{
    for (int i = 0; i < n; i++)
    {
        E1[i] = E[i];

    }
}
/////////////////////////////////////////////////////////Подсчитывает энергию одного спина
float EnergyIJ(int i,short *spin,short *svaz)
{
    float Ener = 0;

    if(i!=0 && i!=1 && i!=n-1 && i!=n-2)
    {
       Ener = -svaz[i] * spin[i] * spin[i+1] - svaz[i-1] * spin[i] * spin[i-1]- svaz[n+i] * spin[i] * spin[i+2] - svaz[(n+i)-2] * spin[i] * spin[i-2] - Polevar;
    }
    else
    {
        if(i==0)
            Ener = -svaz[0] * spin[0] * spin[1] - svaz[n-1] * spin[0] * spin[n-1] - svaz[n] * spin[0] * spin[2] - svaz[2*n-2] * spin[0] * spin[n-2] - Polevar;
        if(i==1)
            Ener = -svaz[0] * spin[0] * spin[1] - svaz[1] * spin[1] * spin[2] - svaz[n+1] * spin[1] * spin[3] - svaz[2*n-1] * spin[1] * spin[n-1] -  Polevar;
        if(i==n-1)
            Ener = -svaz[n-1] * spin[n-1] * spin[0] - svaz[n-2] * spin[n-1] * spin[n-2] - svaz[2*n-1] * spin[n-1] * spin[1] - svaz[n-1] * spin[n-1] * spin[0] - Polevar;
        if(i==n-2)
            Ener = -svaz[n-2] * spin[n-2] * spin[n-1] - svaz[n-3] * spin[n-3] * spin[n-2] - svaz[2*n-2] * spin[n-2] * spin[0] - svaz[2*n-4] * spin[n-2] * spin[n-4] - Polevar;
    }



    return Ener;
}
void MCProhod(double temp,short *spin,short *svaz,float *E)
{

    float En1,En2;
    int step=0;
    double slch;
    double veroyatnost=0;
    while(step<=(n))
    {
        int slspin=rand()%(n);
        int i = slspin;
        int var = spin[i];
        En1 = EnergyIJ(i, spin, svaz);
        spin[i] *= -1;
        En2 = EnergyIJ(i, spin, svaz);
        slch=((double)rand())/RAND_MAX;
        if (En2 <= En1)
            veroyatnost = 1;
        else
            veroyatnost = exp(-(En2 - En1) / temp);
        if (slch < veroyatnost)
            E[i] = En2;
        else
            spin[i] *= -1;
        step++;
    }
    Energy(spin,svaz,E);
}


///////////////////////////////////////////////////////
void MCProhod2(double temp,short *spin,short *svaz,float *E, int rank)
{

    float En1,En2;
    int step=0;
    double slch;
    double veroyatnost=0;

        int slspin=rand()%(n);
        int i = slspin;
        int var = spin[i];
        En1 = EnergyIJ(i, spin, svaz);
        spin[i] *= -1;
        En2 = EnergyIJ(i, spin, svaz);
        slch=((double)rand())/RAND_MAX;
//        slch=myrandom(my_random_engine);
        if (En2 <= En1)
            veroyatnost = 1;
        else
            veroyatnost = exp(-(En2 - En1) / temp);
        if (slch < veroyatnost)
        {
            E[i] = En2;
                if (i == 0)
                {
                   E[n-2]= EnergyIJ(n-2,spin,svaz);
                   E[n-1]= EnergyIJ(n-1,spin,svaz);
                   E[1]= EnergyIJ(1,spin,svaz);
                   E[2]= EnergyIJ(2,spin,svaz);
                }
                if (i == 1)
                {
                   E[n-1]= EnergyIJ(n-1,spin,svaz);
                   E[0]= EnergyIJ(0,spin,svaz);
                   E[2]= EnergyIJ(2,spin,svaz);
                   E[3]= EnergyIJ(3,spin,svaz);

                }

                if (i == n - 1)
                {
                    E[1] = EnergyIJ(1,spin,svaz);
                    E[0] = EnergyIJ(0,spin,svaz);
                    E[n-2] = EnergyIJ(n-2,spin,svaz);
                    E[n-3] = EnergyIJ(n-3,spin,svaz);
                }

                if (i == n - 2)
                {
                    E[0] = EnergyIJ(0,spin,svaz);
                    E[n-1] = EnergyIJ(n-1,spin,svaz);
                    E[n-3] = EnergyIJ(n-3,spin,svaz);
                    E[n-4] = EnergyIJ(n-4,spin,svaz);
                }

                if (i != 0 && i != 1 && i != n-2 && i != n - 1)
                {
                        E[i-2] = EnergyIJ(i-2,spin,svaz);
                        E[i-1] = EnergyIJ(i-1,spin,svaz);
                        E[i+1] = EnergyIJ(i+1,spin,svaz);
                        E[i+2] = EnergyIJ(i+2,spin,svaz);
                }


        }
        else
            spin[i] *= -1;
        step++;


//      Energy(spin,svaz,E);


}
///////////////////////////////////////////////////////////Подсчет намагниченности
double MagnMet(short *spin)
{
    double Magnvar=0;
    for (int i = 0; i < n; ++i)
    {
        Magnvar += spin[i];

    }
    return abs(Magnvar/(n));
}
////////////////////////////////////////////////Для подсчета макс.кластера для спинового стекла
int Claster(float *E1)
{
   int tmp = 0;
    int t;
    int Max = 0;
    int *Ochered;
    Ochered = new int[n];
    for (int i = 0; i < n; i++)
    {
        if (E1[i] == -4 || E1[i] == -2)
            {
                tmp = 1;
                int top = i;
                t = MaxClass(top, E1, tmp,Ochered);
                if (t >= Max)
                    Max = t;
                t = 0;
            }

    }
    delete [] Ochered;
    Ochered = NULL;
    return Max;
}
//////////////////////////////////////////////////////Для подсчета макс.кластера для ферромагеника
int ClasterFerr(float *E1)
{
    int tmp = 0;
    int t;
    int Max = 0;
    int *Ochered;
    Ochered = new int[n];
    for (int i = 0; i < n; i++)
    {
        if (E1[i] == -4 )
            {
                tmp = 1;
                int top = i;
                t = MaxClassFerr(top, E1, tmp, Ochered);
                if (t >= Max)
                    Max = t;
                t = 0;
            }

    }
    delete [] Ochered;
    Ochered = NULL;
    return Max;
}
//////////////////////////////////////////////////////Для подсчета макс.кластера для спинового стекла
int MaxClass(int per,float *E1,int tmp,int *Ochered)
{
    int sl, sp, slsl, spsp, top;
    int w = 0;
    int r = 1;
    Ochered[0] = per;
    E1[per] = 33;
    while (w < r)
    {
        top = Ochered[w];

        ////
        sl = sosed_sl(top);
        if (E1[sl] == -4 || E1[sl] == -2)
        {
            Ochered[r] = sl;
            E1[sl] = 33;
            r++;
        }
        ////
        sp = sosed_sp(top);
        if (E1[sp] == -4 || E1[sp] == -2)
        {
            Ochered[r] = sp;
            E1[sp] = 33;
            r++;
        }
        spsp = sosed_spsp(top);
        if (E1[spsp] == -4 || E1[spsp] == -2)
        {
            Ochered[r] = spsp;
            E1[spsp] = 33;
            r++;
        }
        slsl = sosed_slsl(top);
        if (E1[slsl] == -4 || E1[slsl] == -2)
        {
            Ochered[r] = slsl;
            E1[slsl] = 33;
            r++;
        }
        ////
        w++;

    }
    return r;


}
/////////////////////////////////////////////////////Для подсчета макс.кластера для ферромагеника
int MaxClassFerr(int per,float *E1,int tmp, int *Ochered)
{
     int sl, sp, slsl, spsp,top;
    int w = 0;
    int r = 1;
    Ochered[0] = per;
    E1[per] = 33;
    while (w < r)
    {
        top = Ochered[w];

        ////
        sl = sosed_sl(top);
        if (E1[sl] == -4 )
        {
            Ochered[r] = sl;
            E1[sl] = 33;
            r++;
        }
        ////
        sp = sosed_sp(top);
        if (E1[sp] == -4 )
        {
            Ochered[r] = sp;
            E1[sp] = 33;
            r++;
        }
        spsp = sosed_spsp(top);
        if (E1[spsp] == -4 )
        {
            Ochered[r] = spsp;
            E1[spsp] = 33;
            r++;
        }
        slsl = sosed_slsl(top);
        if (E1[slsl] == -4 )
        {
            Ochered[r] = slsl;
            E1[slsl] = 33;
            r++;
        }
        ////
        w++;

    }
    return r;
}
//////////////////////////////////////////////////// Служебные функции для посчета мах.кластера

int sosed_sl( int top)
{
    int i = top;
    int sl = 0;
    if (i == 0)
        sl = n-1;
    else
        sl = i - 1;
    return sl;
}
int sosed_sp( int top)
{
    int i = top;
    int sp = 0;
    if (i == n - 1)
        sp = 0;
    else
        sp = i + 1;
    return sp;
}
int sosed_slsl( int top)
{
    int i = top;
    int slsl = 0;
    if (i == 0)
       slsl=n-2;
    else if (i == 1)
       slsl=n-1;
    else
       slsl = i - 2;
    return slsl;
}
int sosed_spsp( int top)
{
    int i = top;
    int spsp = 0;
    if (i == n - 2)
            spsp = 0;
    else if (i == n - 1)
            spsp = 1;
    else
        spsp = i + 2;

    return spsp;
}
//////////////////////////////////////////////////////////Суммирует энергию системы
float SumEnergy( float *E)
{
    float Sum = 0;
    for (int i = 0; i < n; i++)
    {
        Sum = Sum + E[i];
    }
    return Sum/2;
}
/////////////////////////////////////////////////////////// Получает поле
void POLE(short *spin,short *svaz,float *E1)
{
    int x = n;
    E1[0] = -svaz[0] * spin[1] - svaz[n-1] * spin[n-1] - Polevar;
    E1[n-1] = -svaz[n-2] * spin[n-2] - svaz[n-1] * spin[0] - Polevar;
    for (int i = 1; i < n - 1; i++)
    {
        E1[i] = -svaz[i] * spin[i+1] - svaz[i-1] * spin[i-1] - Polevar;
    }

}
///////////////////////////////////////////////////////////////
void PhazDiagramm(int *var02,int *var00,int *var20, float *E)
{
    *var02=0;
    *var00=0;
    *var20=0;
    for (int i = 0; i < n; i++)
    {
        if(E[i]==-2)
        {
            *var02+=1;
        } else if (E[i]==0)
        {
            *var00+=1;
        } else if (E[i]==2)
        {
            *var02+=1;
        }

    }

}

static void SROTKL(int size,double Outputmass[],double *SRZnach,double *SRKVOtkl)
{
    *SRZnach = 0;
    *SRKVOtkl = 0;
    for (int i = 0; i < size; i++)
    {
        *SRZnach += Outputmass[i];
    }
    *SRZnach = *SRZnach / size;
    for (int i = 0; i < size; i++)
    {
        Outputmass[i] = pow((Outputmass[i] - *SRZnach), 2);
        *SRKVOtkl += Outputmass[i];
    }
    *SRKVOtkl = sqrt(*SRKVOtkl / size);

}

static void SROTKL(int size,int Outputmass[],double *SRZnach,double *SRKVOtkl)
{
    *SRZnach = 0;
    *SRKVOtkl = 0;
    for (int i = 0; i < size; i++)
    {
        *SRZnach += Outputmass[i];
    }
    *SRZnach = *SRZnach / size;
    for (int i = 0; i < size; i++)
    {
        Outputmass[i] = pow((Outputmass[i] - *SRZnach), 2);
        *SRKVOtkl += Outputmass[i];
    }
    *SRKVOtkl = sqrt(*SRKVOtkl / size);

}

int main(int argc, char **argv)
{    
    int rank;
    int size;
        cout<<"START:"<<endl;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
srand((unsigned)time(NULL)+rank);
 int CO = Prohod_MC_sampling/CountOut;
   // int CO = prohod_MC/CountOut;

    ///////// Создание массива spin[n]////////////
    short *spin=new short[n];
    /////////////////////////////////////

    //Создание массива svaz[n]////////////
    short *svaz = new short[2*n];
    /////////////////////////////////////////

    ////// Создание массива E[n]////////////

    float *E=new float[n];
    /////////////////////////////////////

    ////// Создание массива E1[n]////////////
    float *E1=new float[n];
    ////////////////////////////////////

    ////////// Задание связей ферромагнетика///////////
   for(int i=0;i<2*n;++i)
        {
            svaz[i]=1;

        }
    /////////////массивы для вывода//////////////////////////

        double *recvAE = new double[size];
    double *recvAE2 = new double[size];
    double *recvAE4 = new double[size];
    double *recvAM = new double[size];
    double *recvAM2 = new double[size];
    double *recvAM4 = new double[size];
    double *recvAPP = new double[size];
    double *recvAPPF = new double[size];

    double *recvHeatCapacity = new double[size];
    double *recvVospr = new double[size];

    double *recvBCenergy = new double[size];
    double *recvBCmagn = new double[size];
    double *recvBCPPF = new double[size];
    double *recvBCPP = new double[size];
    



    //ofstream outE("energy.dat",ios::out);           //энергия
    ofstream outAE("energyAVG.dat",ios::out);       //средняя энергия
    ofstream outAE2("energyAVG2.dat",ios::out);       //средняя энергия квадрат
    ofstream outAE4("energyAVG4.dat",ios::out);       //средняя энергия 4 степень
 
ofstream outHeatCapacity("C.dat",ios::out);     //Теплоемкость
ofstream outVospr("X.dat",ios::out);            //Восприимчивость
    ofstream outAPP("APP.dat",ios::out);            //Средний параметр порядка 1
    ofstream outAPPF("APPF.dat",ios::out);          //Средний параметр порядка 2
    ofstream outAM("AMagn.dat",ios::out);           //Средняя намагниченность
    ofstream outAM2("AMagn2.dat",ios::out);         
    ofstream outAM4("AMagn4.dat",ios::out);
    // биндеры
    ofstream outBCenergy("BC_energy.dat",ios::out); //биндер по энергии
    ofstream outBCmagn("BC_magn.dat",ios::out);     //биндер по намагниченности
    ofstream outBCPP("BC_PP.dat",ios::out);         //биндер по параметру порядка
    ofstream outBCPPF("BC_PPF.dat",ios::out);       //биндер по параметру порядка2

 
    ////////////////////////////
    SpinFerr(spin);
    //srand((unsigned)time(NULL)+rank);
    //Slspin(spin);
    Energy(spin,svaz,E);
    CopyE(E,E1);
    /////////////////// время///////////////
    
    //////MC///////////////////////////////


    double SumenergyVar,MagnVar,MagnVar2,MagnVar4,Teploemkost,Vospriimchivost;
    double predsumenergy=0;
    double predMagn=0;
    double maxSF,PPF,maxS,PP;
    int var02=0;
    int var00=0;
    int var20=0;

    int varPOLE02=0;
    int varPOLE00=0;
    int varPOLE20=0;
    double SumenergyVar2, SumenergyVar4;
    double Aenergy2=0;
    double Aenegry4=0;
    double Amagn2=0;
    double Amagn4=0;
    double APPF2=0;
    double APPF4=0;
    double APP2=0;
    double APP4=0;
    double X_o=0;


    //    double temp = rank*Dt;int Prohod=0;
    //    int cicle=1;

    vector<double> avgValue;    // AE
    vector<double> avgValue2;    // AE2
    vector<double> avgValue4;    // AE4
    vector<double> avgValuePP;  // APP
    vector<double> avgValuePPF; // APPF
    vector<double> avgValueMagn;// AM

    vector<double> CValue;        // E
    vector<double> CValuePP;  // PP
    vector<double> CValuePPF; // PPF
    vector<double> CValueMagn;// M



    double Avalue=0;    // для усреднения энергии
    double Aenergy=0;
    unsigned long long index=0;


    double AvaluePP=0;  // для усреднения PP
    double temPP;
    double AvaluePPF=0; // для усреднения PPF
    double temPPF;
    double AvalueMagn=0; // для усреднения Magn
    double temMagn;

        //++++++++++++++++++++++++++++++++++++++
        double e_AVG=0;
        double PP1_AVG=0;
        double e1=0;
        double pp1=0;

        // для вывода БК for PP2

        double RBCenergy=0;
        double RBCmagn=0;
        double RBCPPF=0;
        double RBCPP=0;
     
        double E_2=0;
        double C_o=0;
       
        //++++++++++++++++++++++++++++++++++++++
        //для репличного обмена
    double exchange_t=0; 
    double exchange_e=0;  
    int yes_no_exchange=0; 

    double probability_of_exchange=0; 
 //   bool replica_on=1;
    MPI_Status status; 
    short *exchange_state=new short[n];
    


    //++++++++++++++++++++++++++++++++++++++

    double temp=(rank+1)*(maxtemp-mintemp)/size;
    long int Prohod=0;
    while(Prohod<Prohod_MC_equilibration)
    {
        e1=SumEnergy(E);
        if(replica_on)
        {   
            if(Prohod%between_exchange_equilibration==0)
            {
                for(int i=size-1;i>0;--i)
                {
                    if(rank==i)
                    {
                        MPI_Send(&temp,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD);
                        MPI_Send(&e1,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD);
                        MPI_Recv(&yes_no_exchange,1, MPI_INT, i-1, 0, MPI_COMM_WORLD,&status);
                        if(!yes_no_exchange);

                        else
                        {
                            MPI_Send(spin,n, MPI_SHORT, i-1, 0, MPI_COMM_WORLD);
                            MPI_Recv(exchange_state,n, MPI_SHORT, i-1, 0, MPI_COMM_WORLD,&status);
                            MPI_Recv(&e1,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD,&status);
                            CopyE(exchange_state,spin);
                        }
                    }
                        if(rank==i-1)
                        {
                            MPI_Recv(&exchange_t,1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,&status);
                            MPI_Recv(&exchange_e,1, MPI_DOUBLE, i, 0,MPI_COMM_WORLD,&status);
                            probability_of_exchange=exp(-((1/temp)-(1/exchange_t))*(exchange_e-e1));
                            if (probability_of_exchange<(double)(rand())/RAND_MAX)
                            {
                                yes_no_exchange=0;
                                MPI_Send(&yes_no_exchange,1, MPI_INT, i, 0, MPI_COMM_WORLD);

                            }
                            else
                            {
                                yes_no_exchange=1;
                                MPI_Send(&yes_no_exchange,1, MPI_INT, i, 0, MPI_COMM_WORLD); //порядок?
                                MPI_Recv(exchange_state,n, MPI_SHORT, i, 0, MPI_COMM_WORLD,&status);
                                MPI_Send(spin,n, MPI_SHORT, i, 0, MPI_COMM_WORLD); 
                                MPI_Send(&e1,1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                                CopyE(exchange_state, spin);
                            }
                        }

                        MPI_Barrier(MPI_COMM_WORLD);
                }
            }
        }

        Energy(spin,svaz,E);
        MCProhod2(temp,spin,svaz,E,rank);    
         if(rank==0)
        {
            if(Prohod%(Prohod_MC_equilibration/100)==0) //было if(Prohod%(Prohod_MC_sampling/100)==0)
                cout<<"Equilibration status: "<<(Prohod)/(Prohod_MC_equilibration/100)+1<<endl; // было cout<<"Status: "<<Prohod/(Prohod_MC_sampling/100)<<endl;
        }
        Prohod++;
    }
    Prohod=0;

 ///////////////////////////////////   
    while(Prohod<Prohod_MC_sampling)
    {
        e1=SumEnergy(E);
        if(replica_on)
        {   
            if(Prohod%between_exchange_equilibration==0)
            {
                for(int i=size-1;i>0;--i)
                {
                    if(rank==i)
                    {
                        MPI_Send(&temp,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD);
                        MPI_Send(&e1,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD);
                        MPI_Recv(&yes_no_exchange,1, MPI_INT, i-1, 0, MPI_COMM_WORLD,&status);
                        if(!yes_no_exchange);

                        else
                        {
                            MPI_Send(spin,n, MPI_SHORT, i-1, 0, MPI_COMM_WORLD);
                            MPI_Recv(exchange_state,n, MPI_SHORT, i-1, 0, MPI_COMM_WORLD,&status);
                            MPI_Recv(&e1,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD,&status);
                            CopyE(exchange_state,spin);
                        }
                    }
                        if(rank==i-1)
                        {
                            MPI_Recv(&exchange_t,1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,&status);
                            MPI_Recv(&exchange_e,1, MPI_DOUBLE, i, 0,MPI_COMM_WORLD,&status);
                            probability_of_exchange=exp(-((1/temp)-(1/exchange_t))*(exchange_e-e1));
                            if (probability_of_exchange<(double)(rand())/RAND_MAX)
                            {
                                yes_no_exchange=0;
                                MPI_Send(&yes_no_exchange,1, MPI_INT, i, 0, MPI_COMM_WORLD);

                            }
                            else
                            {
                                yes_no_exchange=1;
                                MPI_Send(&yes_no_exchange,1, MPI_INT, i, 0, MPI_COMM_WORLD); //порядок?
                                MPI_Recv(exchange_state,n, MPI_SHORT, i, 0, MPI_COMM_WORLD,&status);
                                MPI_Send(spin,n, MPI_SHORT, i, 0, MPI_COMM_WORLD); 
                                MPI_Send(&e1,1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                                e1=exchange_e;
                                CopyE(exchange_state, spin);
                            }
                        }

                        MPI_Barrier(MPI_COMM_WORLD);
                }
            }
        }

        Energy(spin,svaz,E);
        MCProhod2(temp,spin,svaz,E,rank); 
        if(Prohod%CO==0)
        {
            index++;
             if (Polevar==0)
            {
            CopyE(E,E1);                    // при расчете кластера, коцаеца массив энергий, поэтому мы его копируем и работаем с копией
            maxSF = ClasterFerr(E1);
            PPF = maxSF / (n);                  // 
            temPPF=(PPF+(index-1)*temPPF)/(index);      // не забыть помножить !!!*N куммулятив ферромагнитного кластера
            APPF2=(PPF+(index-1)*temPPF)/(index);
            APPF4=(PPF+(index-1)*temPPF)/(index);

            CopyE(E,E1);                        // при расчете кластера, коцаеца массив энергий, поэтому мы его копируем и работаем с копией
            maxS = Claster(E1);
            PP= maxS / (n);
            temPP=(PP+(index-1)*temPP)/(index);      // !!!*N куммулятив смешанного кластера
            APP2=(PP+(index-1)*APP2)/(index);
            APP4=(PP+(index-1)*APP4)/(index);
            }
            SumenergyVar = SumEnergy(E);     /////////////
            SumenergyVar2=pow(SumenergyVar,2);
            SumenergyVar4=pow(SumenergyVar,4);
            Aenergy=(SumenergyVar+(index-1)*Aenergy)/(index);   //было tem213=Avalue*(CO)/(Prohod);  
            Aenergy2=(SumenergyVar2+(index-1)*Aenergy2)/(index);
            Aenegry4=(SumenergyVar4+(index-1)*Aenegry4)/(index);


            MagnVar=MagnMet(spin);  
            MagnVar2=MagnVar*MagnVar;
            MagnVar4=MagnVar*MagnVar*MagnVar*MagnVar;       /////////////
            temMagn=(MagnVar+(index-1)*temMagn)/(index);   // !!!*N куммулятив смешанного кластера
            Amagn2=(MagnVar2+(index-1)*Amagn2)/(index);
            Amagn4=(MagnVar4+(index-1)*Amagn4)/(index);


            ////////////////////БК_И_Теплоемкость////////////////////
            e_AVG+=SumenergyVar;
            E_2+=SumenergyVar*SumenergyVar;
            PP1_AVG+=SumenergyVar; //e_AVG+=SumenergyVar;

            //  }
            if(rank==0)
            {
                if(Prohod%(Prohod_MC_sampling/100)==0) //было if(Prohod%(Prohod_MC_sampling/100)==0)
                    cout<<"Status: "<<(Prohod)/(Prohod_MC_sampling/100)+1<<endl; // было cout<<"Status: "<<Prohod/(Prohod_MC_sampling/100)<<endl;
            }
        }

        Prohod++;
    }
    cout<<"index:"<<index<<"    Coutout:"<<CountOut<<endl;
    C_o=((Aenergy2)-(Aenergy*Aenergy))/(temp*temp*n); //  теплоёмкость //  теплоёмкость
    X_o=((Amagn2)-(temMagn*temMagn))/(temp*n);
     RBCenergy = 1-(Aenegry4/(3*pow(Aenergy2,2)));  // результирующий БК энергии
    RBCmagn = 1-(Amagn4/(3*pow(Amagn2,2)));  // результирующий БК магнитный
    RBCPPF = 1-(APPF4/(3*pow(APPF2,2)));  // результирующий БК параметра порядка 2
    RBCPP = 1-(APP4/(3*pow(APP2,2)));  // результирующий БК параметра порядка 


    MPI_Gather(&X_o, 1, MPI_DOUBLE, recvVospr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&C_o, 1, MPI_DOUBLE, recvHeatCapacity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Теплоемкость
    MPI_Gather(&RBCenergy, 1, MPI_DOUBLE, recvBCenergy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК энергию
    MPI_Gather(&RBCmagn, 1, MPI_DOUBLE, recvBCmagn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК магн
    MPI_Gather(&RBCPPF, 1, MPI_DOUBLE, recvBCPPF, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК PPF
    MPI_Gather(&RBCPP, 1, MPI_DOUBLE, recvBCPP, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК PP
   MPI_Gather(&Aenergy, 1, MPI_DOUBLE, recvAE, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&Aenergy2, 1, MPI_DOUBLE, recvAE2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&Aenegry4, 1, MPI_DOUBLE, recvAE4, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&temMagn, 1, MPI_DOUBLE, recvAM, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     MPI_Gather(&Amagn2, 1, MPI_DOUBLE, recvAM2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&Amagn4, 1, MPI_DOUBLE, recvAM4, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&temPP, 1, MPI_DOUBLE, recvAPP, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&temPPF, 1, MPI_DOUBLE, recvAPPF, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     if(rank==0) //countOut - количество выводов!!! CO- какое значение из общего количества
    {
        for(int i=0; i<size; ++i)
        {
            for(int j=0; j<CountOut; ++j)
            {

                // outE<<(i+1)*((maxtemp-mintemp)/size)<<'\t'<<recvE[i*CountOut+j]<<endl;
                
                // outPP<<(i+1)*((maxtemp-mintemp)/size)<<'\t'<<recvPP[i*CountOut+j]<<endl;
                // outPPF<<(i+1)*((maxtemp-mintemp)/size)<<'\t'<<recvPPF[i*CountOut+j]<<endl;                
            }
             outAPPF<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAPPF[i]<<endl;
            outAPP<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAPP[i]<<endl;
            outAM<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAM[i]<<endl;
             outAM2<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAM2[i]<<endl;
              outAM4<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAM4[i]<<endl;
            outAE<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAE[i]<<endl;
            outAE2<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAE2[i]<<endl;
            outAE4<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAE4[i]<<endl;

            outHeatCapacity<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvHeatCapacity[i]<<endl;                 // Теплоемкость

            outVospr<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvVospr[i]<<endl;  //воспириимчивость
            outBCenergy<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCenergy[i]<<endl; // БК energy
            outBCmagn<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCmagn[i]<<endl; // БК magn
            outBCPPF<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCPPF[i]<<endl; // БК PPF
            outBCPP<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCPP[i]<<endl; // БК PP
       }
    }
    
    delete []svaz;
    delete []E;
    delete []E1;
    delete []spin;
    delete []exchange_state;

    delete []recvAPPF;
    delete []recvAPP;
    delete []recvAM;
  
    delete []recvAM2;
    delete []recvAM4;
    delete []recvAE;
    delete []recvAE2;
    delete []recvAE4;
    delete []recvBCenergy;
    delete []recvBCmagn;
    delete []recvBCPPF;
    delete []recvBCPP;
    delete []recvHeatCapacity;

    MPI_Finalize();
    cout<<"FINISH"<<endl;

    return 0;
}