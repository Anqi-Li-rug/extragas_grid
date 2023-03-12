/*
 *  fitsutil.cpp
 *  Bulge
 *
 *  Created by Antonino Marasco on 16/04/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#define PI 3.1415926535

#include "fitsutil.h"
#include "cmath"
#include "zap.h"

extern bool verbose;

inline double gauss1D(double x, double sigma)
{
    return exp(-x*x/(2.*sigma*sigma))/(sigma*sqrt(2.*PI));
};

inline double gauss2D(double x, double y, double sigmax, double sigmay)
{
    return gauss1D(x,sigmax)*gauss1D(y,sigmay);
};

inline double gauss3D(double x, double y, double z, double sigmax, double sigmay, double sigmaz)
{
    return gauss1D(x,sigmax)*gauss1D(y,sigmay)*gauss1D(z,sigmaz);
};

////////////////////////////
//         MAPS           //
////////////////////////////
void fitsutil::read_map(char* input_fits, vector<double> &x, vector<double> &y, double* &map_ptr)
{
    //streambuf* orig_buf = cout.rdbuf();
    //cout.rdbuf(NULL);
    fitsfile *fptr; //puntatore a oggetto fitsfile
    int status = 0; //inizializzo lo stato di fptr
    int nfound, i, j;
    long int naxes[2]; //array di 2 elementi contenenti il numero di dati per ogni dimensione
    double x_value, x_delta, y_value, y_delta, x_pix=1, y_pix=1;
    //apertura del file
    if (fits_open_file(&fptr, input_fits, READONLY, &status))
    {
        cout<<"ERROR IN CFITSIO!"<<endl;
        fits_report_error(stderr, status);
    }
    //lettura della keyword NAXIS del file namefile
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status))  fits_report_error(stderr, status);
    // inizializzazione delle dimensioni della matrice di dati
    x.resize(naxes[0]);
    y.resize(naxes[1]);
    long int dim = x.size()*y.size();
    // lettura delle keyword contenenti i valori iniziali dei vettori e il loro incremento. si suppone siano le stesse per entrambi i files
    printf("Reading %i x %i fits file \n", naxes[0], naxes[1]);
    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &x_value, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &y_value, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TDOUBLE, "CDELT1", &x_delta, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TDOUBLE, "CDELT2", &y_delta, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TDOUBLE, "CRPIX1", &x_pix, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TDOUBLE, "CRPIX2", &y_pix, NULL, &status))  fits_report_error(stderr, status);
    //chiusura dei files
    if (fits_close_file(fptr, &status))  fits_report_error(stderr, status);
    // riempimento dei vettori x e y
    for(i=0;i<x.size();i++) x[i]=x_value+(i+1-x_pix)*x_delta;
    for(j=0;j<y.size();j++) y[j]=y_value+(j+1-y_pix)*y_delta;
    
    zaparr(map_ptr);
    double nullval = 0;
    int anynull;
    if (fits_open_file(&fptr, input_fits, READONLY, &status))  fits_report_error(stderr, status); //apertura file
    map_ptr = new double[dim];
    if (fits_read_img(fptr, TDOUBLE, 1, dim, &nullval, map_ptr, &anynull, &status))  fits_report_error(stderr, status);
    if (fits_close_file(fptr, &status))  fits_report_error(stderr, status); //chiusura file
}

void fitsutil::print_map(char* output_fits, vector<double> &x, vector<double> &y, double* &map_ptr, char* ctype1, char* ctype2, char* cunit1, char* cunit2, char* bunit)
{
    FILE* fitsheader=NULL;
    fitsheader=fopen("myheader.txt", "w");
    fprintf(fitsheader,"CRVAL1 = %f \n",x[0]);
    fprintf(fitsheader,"CRVAL2 = %f \n",y[0]);
    fprintf(fitsheader,"CDELT1 = %f \n",x[1]-x[0]);
    fprintf(fitsheader,"CDELT2 = %f \n",y[1]-y[0]);
    fprintf(fitsheader,"CTYPE1 = %s \n",ctype1);
    fprintf(fitsheader,"CTYPE2 = %s \n",ctype2);
    fprintf(fitsheader,"CUNIT1 = %s \n",cunit1);
    fprintf(fitsheader,"CUNIT2 = %s \n",cunit2);
    fprintf(fitsheader,"CRPIX1 = %i \n",1);
    fprintf(fitsheader,"CRPIX2 = %i \n",1);
    fprintf(fitsheader,"BUNIT  = %s \n",bunit);
    char btype[]="intensity";
    fprintf(fitsheader,"BTYPE  = %s \n",btype);
    fclose(fitsheader);
    fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    long  fpixel = 1, naxis = 2, nelements;
    long dnaxes[2] = {static_cast<long>(x.size()), static_cast<long>(y.size())};
    long naxes[2] = {(int) y.size(), (int) x.size()};
    remove(output_fits);               /* Delete old file if it already exists */
    int status = 0;    /* initialize status before calling fitsio routines */
    fits_create_file(&fptr, output_fits, &status);   /* create new file */
    /* Create the primary array image (32-bit floating pixels */
    fits_create_img(fptr, FLOAT_IMG, naxis, dnaxes, &status);
    fitsheader=fopen("myheader.txt", "r");
    char key[20]; //7
    char eq[20]; //3
    char cval[20]; //10
    char c = 'a';
    while (c!=EOF)
    {
        fscanf(fitsheader, "%s %s %s",key,eq,cval);
        float val=atof(cval);
        c=fgetc(fitsheader);
        if (val == 0) fits_update_key(fptr, TSTRING, key, &cval, " ", &status);
        if (val != 0) fits_update_key(fptr, TFLOAT, key, &val, " ", &status);
    }
    fclose(fitsheader);
    
    printf("Writing %i x %i fits file \n", naxes[1], naxes[0]);
    nelements = naxes[0] * naxes[1]; /* number of pixels to write */
    fits_write_img(fptr, TDOUBLE, fpixel, nelements, map_ptr, &status);
    fits_close_file(fptr, &status);            /* close the file */
    fits_report_error(stderr, status);  /* print out any error messages */
    cout<<"File "<<output_fits<<" written."<<endl;
}


////////////////////////
//       CUBES        //
////////////////////////

void fitsutil::read_cube(char* input_fits, vector<double> &x, vector<double> &y, vector<double> &z, int &crpix1, int &crpix2, int &crpix3, double* &cube_ptr)
{
    //streambuf* orig_buf = cout.rdbuf();
    //cout.rdbuf(NULL);
    fitsfile *fptr; //puntatore a oggetto fitsfile
    int status = 0; //inizializzo lo stato di fptr
    int nfound, i, j, k;
    long int naxes[3]; //array di 3 elementi contenenti il numero di dati per ogni dimensione
    double x_value, x_delta, y_value, y_delta, z_value, z_delta;
    //apertura del file
    if (fits_open_file(&fptr, input_fits, READONLY, &status))
    {
        cout<<"ERROR IN CFITSIO!"<<endl;
        fits_report_error(stderr, status);
    }
    //lettura della keyword NAXIS del file namefile
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 3, naxes, &nfound, &status))  fits_report_error(stderr, status);
    // inizializzazione delle dimensioni della matrice di dati
    x.resize(naxes[0]);
    y.resize(naxes[1]);
    z.resize(naxes[2]);
    long int dim = x.size()*y.size()*z.size();
    
    // lettura delle keyword contenenti i valori iniziali dei vettori e il loro incremento. si suppone siano le stesse per entrambi i files
    //printf("Reading %i x %i x %i fits file \n", naxes[0], naxes[1], naxes[2]);
    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &x_value, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &y_value, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TDOUBLE, "CRVAL3", &z_value, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TDOUBLE, "CDELT1", &x_delta, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TDOUBLE, "CDELT2", &y_delta, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TDOUBLE, "CDELT3", &z_delta, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TINT, "CRPIX1", &crpix1, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TINT, "CRPIX2", &crpix2, NULL, &status))  fits_report_error(stderr, status);
    if (fits_read_key(fptr, TINT, "CRPIX3", &crpix3, NULL, &status))  fits_report_error(stderr, status);
    
    //chiusura dei files
    if (fits_close_file(fptr, &status))  fits_report_error(stderr, status);
    // riempimento dei vettori x e y
    for(i=0;i<x.size();i++) x[i]=x_value + (i+1-crpix1)*x_delta;
    for(j=0;j<y.size();j++) y[j]=y_value + (j+1-crpix2)*y_delta;
    for(k=0;k<z.size();k++) z[k]=z_value + (k+1-crpix3)*z_delta;
    
    //zaparr(cube_ptr);
    double nullval = NAN;
    int anynull;
    if (fits_open_file(&fptr, input_fits, READONLY, &status))  fits_report_error(stderr, status); //apertura file
    cube_ptr = new double [dim];
    if (fits_read_img(fptr, TDOUBLE, 1, dim, &nullval, cube_ptr, &anynull, &status))  fits_report_error(stderr, status);
    if (fits_close_file(fptr, &status))  fits_report_error(stderr, status); //chiusura file
}

void fitsutil::print_cube(char* output_fits, vector<double> &x, vector<double> &y, vector<double> &z, int crpix1, int crpix2, int crpix3, double* &cube_ptr, char* ctype1, char* ctype2, char* ctype3, char* cunit1, char* cunit2, char* cunit3, char* bunit)
{
    //streambuf* orig_buf = cout.rdbuf();
    //cout.rdbuf(NULL);
    FILE* fitsheader=NULL;
    fitsheader=fopen("myheader.txt", "w");
    fprintf(fitsheader,"CRVAL1 = %f \n",x[crpix1-1]);
    fprintf(fitsheader,"CRVAL2 = %f \n",y[crpix2-1]);
    fprintf(fitsheader,"CRVAL3 = %f \n",z[crpix3-1]);
    fprintf(fitsheader,"CDELT1 = %f \n",x[1]-x[0]);
    fprintf(fitsheader,"CDELT2 = %f \n",y[1]-y[0]);
    fprintf(fitsheader,"CDELT3 = %f \n",z[1]-z[0]);
    fprintf(fitsheader,"CTYPE1 = %s \n",ctype1);
    fprintf(fitsheader,"CTYPE2 = %s \n",ctype2);
    fprintf(fitsheader,"CTYPE3 = %s \n",ctype3);
    fprintf(fitsheader,"CUNIT1 = %s \n",cunit1);
    fprintf(fitsheader,"CUNIT2 = %s \n",cunit2);
    fprintf(fitsheader,"CUNIT3 = %s \n",cunit3);
    fprintf(fitsheader,"CRPIX1 = %i \n",crpix1);
    fprintf(fitsheader,"CRPIX2 = %i \n",crpix2);
    fprintf(fitsheader,"CRPIX3 = %i \n",crpix3);
    fprintf(fitsheader,"EPOCH  = %i \n",2000);
    fprintf(fitsheader,"BUNIT  = %s \n",bunit);
    char btype[]="intensity";
    fprintf(fitsheader,"BTYPE  = %s \n",btype);
    fclose(fitsheader);
    fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    long  fpixel = 1, naxis = 3, nelements;
    /*  long naxes[3] = {RAsize, DECsize, VELsize};
     float array[RAsize][DECsize][VELsize]; */
    long dnaxes[3] = {(int) x.size(), (int) y.size(), (int) z.size()};
    long naxes[3] = {(int) z.size(), (int) y.size(), (int) x.size()};
    remove(output_fits);               /* Delete old file if it already exists */
    int status = 0;    /* initialize status before calling fitsio routines */
    fits_create_file(&fptr, output_fits, &status);   /* create new file */
    /* Create the primary array image (32-bit floating pixels */
    fits_create_img(fptr, FLOAT_IMG, naxis, dnaxes, &status);
    fitsheader=fopen("myheader.txt", "r");
    char key[20]; //7
    char eq[20]; //3
    char cval[20]; //10
    char c='a';
    while (c!=EOF)
    {
        fscanf(fitsheader, "%s %s %s",key,eq,cval);
        float val=atof(cval);
        c=fgetc(fitsheader);
        if (val == 0) fits_update_key(fptr, TSTRING, key, &cval, " ", &status);
        if (val != 0) fits_update_key(fptr, TFLOAT, key, &val, " ", &status);
    }
    fclose(fitsheader);
    //printf("Writing %i x %i x %i fits file \n", dnaxes[0], dnaxes[1], dnaxes[2]);
    nelements = naxes[0] * naxes[1] * naxes[2]; /* number of pixels to write */
    fits_write_img(fptr, TDOUBLE, fpixel, nelements, cube_ptr, &status);
    fits_close_file(fptr, &status);            /* close the file */
    fits_report_error(stderr, status);  /* print out any error messages */
    cout<<"File "<<output_fits<<" written."<<endl;
}

void fitsutil::smooth_cube(myarray::array3d<double> &cube, vector<double> &x, vector<double> &y, vector<double> &z, double oldx, double oldy, double newx, double newy, double sigma_accuracy)
{
    if(verbose) cout<<"Smoothing algorithm in progress"<<endl;
    if(newx>oldx && newy>oldy)
    {
    //evaluating the convolution beam
        double a2 = 0.5*oldx; double b2 = 0.5*oldy;
        double a0 = 0.5*newx; double b0 = 0.5*newy;
        double D2 = a2*a2 - b2*b2; double D0 = a0*a0 - b0*b0;
        double D1 = sqrt(D0*D0 + D2*D2 - 2*D0*D2);
        double argument,a1,b1;
        argument = 0.5*(a0*a0 + b0*b0 - a2*a2 - b2*b2 + D1);
        if(argument>0) {a1 = sqrt(argument);}
        else
        {
            cout<<"ERROR: The new beam is not correct!"<<endl;
            return;
        }
        argument = 0.5*(a0*a0 + b0*b0 - a2*a2 - b2*b2 - D1);
        if(argument>0) {b1 = sqrt(argument);}
        else
        {
            cout<<"ERROR: The new beam is not correct!"<<endl;
            return;
        }
        double beamx = 2.*a1;
        double beamy = 2.*b1;

        if(verbose) cout<<"Old beam (in arcsec): "<<oldx<<","<<oldy<<endl;
        if(verbose) cout<<"New beam (in arcsec): "<<newx<<","<<newy<<endl;
        if(verbose) cout<<"Convolution beam (in arcsec): "<<beamx<<","<<beamy<<endl;

    //setting various stuff
        const double myconst = 2.354820045;
        const double deg_to_arcsec = 3600.;
        const double sigmax = beamx/deg_to_arcsec/myconst;
        const double sigmay = beamy/deg_to_arcsec/myconst;
        
        int i,j,k;
        const int X_SIZE = x.size();
        const int Y_SIZE = y.size();
        const int Z_SIZE = z.size();
        const double dx = x[1]-x[0];
        const double dy = y[1]-y[0];
        myarray::array3d<double> cubetemp(X_SIZE,Y_SIZE,Z_SIZE);
        
        int i1, i2, j1, j2;
        double new_percent = 0, percent = 0;
        double weight,value,temp;
        int state = 0, newstate;
        for(k=0;k<Z_SIZE;k++)
        {
            newstate = static_cast<int>((k+1)*100/Z_SIZE);
            if (newstate != state)
            {
                //if(verbose) cout<<"Progress : "<<newstate<<" percent"<<endl;
                state = newstate;
            }
    //smoothing in the x-direction
            for(j=0;j<Y_SIZE;j++) for(i=0;i<X_SIZE;i++)
            {
                if(isnan(cube(i,j,k))) {cubetemp(i,j,k)=NAN;}
                else
                {
                    i1 = i - fabs(sigma_accuracy*sigmax/dx) - 0.5;
                    i2 = i + fabs(sigma_accuracy*sigmax/dx) + 0.5;
                    if(i1<0) i1 = 0;
                    if(i2>X_SIZE-1) i2=X_SIZE-1;
                    weight = 0; value = 0;
                    for(int n=i1;n<=i2;n++)
                    {
                        if(!isnan(cube(n,j,k)))
                        {
                            temp = gauss1D(x[n]-x[i],sigmax)*fabs(dx);
                            weight+=temp;
                            value +=cube(n,j,k)*temp;
                        }
                    }
                    cubetemp(i,j,k)=value/weight;
                }
            }
    //smoothing in the y-direction
            for(i=0;i<X_SIZE;i++) for(j=0;j<Y_SIZE;j++)
            {
                if(isnan(cubetemp(i,j,k))) {cube(i,j,k)=NAN;}
                else
                {
                    j1 = j - fabs(sigma_accuracy*sigmay/dy) - 0.5;
                    j2 = j + fabs(sigma_accuracy*sigmay/dy) + 0.5;
                    if(j1<0) j1 = 0;
                    if(j2>Y_SIZE-1) j2=Y_SIZE-1;
                    weight = 0; value = 0;
                    for(int n=j1;n<=j2;n++)
                    {
                        if(!isnan(cubetemp(i,n,k)))
                        {
                            temp = gauss1D(y[n]-y[j],sigmay)*fabs(dy);
                            weight+=temp;
                            value +=cubetemp(i,n,k)*temp;
                        }
                    }
                    cube(i,j,k)=value/weight;
                }
            }
        }
    }
}


void fitsutil::smooth_allsky_cube(myarray::array3d<double> &cube, vector<double> &x, vector<double> &y, vector<double> &z, double oldx, double oldy, double newx, double newy, double sigma_accuracy)
{
    if(verbose) cout<<"Smoothing algorithm in progress"<<endl;
    if(newx>oldx && newy>oldy)
    {
        //evaluating the convolution beam
        double a2 = 0.5*oldx; double b2 = 0.5*oldy;
        double a0 = 0.5*newx; double b0 = 0.5*newy;
        double D2 = a2*a2 - b2*b2; double D0 = a0*a0 - b0*b0;
        double D1 = sqrt(D0*D0 + D2*D2 - 2*D0*D2);
        double argument,a1,b1;
        argument = 0.5*(a0*a0 + b0*b0 - a2*a2 - b2*b2 + D1);
        if(argument>0) {a1 = sqrt(argument);}
        else
        {
            cout<<"ERROR: The new beam is not correct!"<<endl;
            return;
        }
        argument = 0.5*(a0*a0 + b0*b0 - a2*a2 - b2*b2 - D1);
        if(argument>0) {b1 = sqrt(argument);}
        else
        {
            cout<<"ERROR: The new beam is not correct!"<<endl;
            return;
        }
        double beamx = 2.*a1;
        double beamy = 2.*b1;
        
        if(verbose) cout<<"Old beam (in degree): "<<oldx<<","<<oldy<<endl;
        if(verbose) cout<<"New beam (in degree): "<<newx<<","<<newy<<endl;
        if(verbose) cout<<"Convolution beam (in degree): "<<beamx<<","<<beamy<<endl;
        
        //setting various stuff
        const double myconst = 2.354820045;
        const double sigmax_midplane = beamx/myconst;
        const double sigmay = beamy/myconst;
        
        int i,j,k;
        const int X_SIZE = x.size();
        const int Y_SIZE = y.size();
        const int Z_SIZE = z.size();
        const double dx = x[1]-x[0];
        const double dy = y[1]-y[0];
        myarray::array3d<double> cubetemp(X_SIZE,Y_SIZE,Z_SIZE);
        
        int i1, i2, j1, j2;
        double new_percent = 0, percent = 0;
        double weight,value,temp,sigmax;
        int state = 0, newstate;
        for(k=0;k<Z_SIZE;k++)
        {
            newstate = static_cast<int>((k+1)*100/Z_SIZE);
            if (newstate != state)
            {
                //if(verbose)cout<<"Progress : "<<newstate<<" percent"<<endl;
                state = newstate;
            }
            //smoothing in the x-direction
            for(j=0;j<Y_SIZE;j++)
            {
                sigmax = sigmax_midplane/cos(y[j]*0.0174532925);
                //sigmax cannot be arbitrarly large, so...
                if(sigma_accuracy*sigmax>90) sigmax=90;
                for(i=0;i<X_SIZE;i++)
                {
                    if(isnan(cube(i,j,k))) {cubetemp(i,j,k)=NAN;}
                    else
                    {
                        i1 = i - fabs(sigma_accuracy*sigmax/dx) - 0.5;
                        i2 = i + fabs(sigma_accuracy*sigmax/dx) + 0.5;
                        weight = 0; value = 0;
                        for(int n=i1;n<=i2;n++)
                        {
                            if(!isnan(cube(n,j,k)))
                            {
                                temp = gauss1D(x[n]-x[i],sigmax)*fabs(dx);
                                weight+=temp;
                                value +=cube(n,j,k)*temp;
                            }
                        }
                        cubetemp(i,j,k)=value/weight;
                    }
                }
            }
            //smoothing in the y-direction
            for(i=0;i<X_SIZE;i++) for(j=0;j<Y_SIZE;j++)
            {
                if(isnan(cubetemp(i,j,k))) {cube(i,j,k)=NAN;}
                else
                {
                    j1 = j - fabs(sigma_accuracy*sigmay/dy) - 0.5;
                    j2 = j + fabs(sigma_accuracy*sigmay/dy) + 0.5;
                    if(j1<0) j1 = 0; if(j2>=y.size()) j2=y.size()-1;
                    weight = 0; value = 0;
                    for(int n=j1;n<=j2;n++)
                    {
                        if(!isnan(cubetemp(i,n,k)))
                        {
                            temp = gauss1D(y[n]-y[j],sigmay)*fabs(dy);
                            weight+=temp;
                            value +=cubetemp(i,n,k)*temp;
                        }
                    }
                    cube(i,j,k)=value/weight;
                }
            }
        }
    }
}


void fitsutil::smooth3D(myarray::array3d<double> &cube, vector<double> &x, vector<double> &y, vector<double> &z, double beamx, double beamy, double beamz, double sigma_accuracy)
{
    if(verbose) if(verbose) cout<<"Smoothing algorithm in progress"<<endl;
    const double myconst = 2.354820045;
    const double sigmax = beamx/myconst;
    const double sigmay = beamy/myconst;
    const double sigmaz = beamz/myconst;
    
    const int X_SIZE = x.size();
    const int Y_SIZE = y.size();
    const int Z_SIZE = z.size();
    const double dx = x[1]-x[0];
    const double dy = y[1]-y[0];
    const double dz = z[1]-z[0];
    
    if(verbose)
    {
        cout<<"Old beam: ("<<fabs(dx)<<","<<fabs(dy)<<","<<fabs(dz)<<")"<<endl;
        cout<<"New beam: ("<<beamx<<","<<beamy<<","<<beamz<<")"<<endl;
    }
    myarray::array3d<double> cubetemp(X_SIZE,Y_SIZE,Z_SIZE);
    
    int i1, i2, j1, j2, k1, k2;
    double x0, y0, z0;
    double weight,weight_here,value;
    long int counter = 0;
    long int totdim = X_SIZE*Y_SIZE*Z_SIZE;
    int state = 0, newstate;
    
    for(int ii=0;ii<X_SIZE;ii++) for(int jj=0;jj<Y_SIZE;jj++) for(int kk=0;kk<Z_SIZE;kk++)
    {
        newstate = 100*counter/totdim;
        if (newstate != state)
        {
            //if(verbose) cout<<"Progress : "<<newstate<<" percent"<<endl;
            state = newstate;
        }
        
        i1 = ii - fabs(sigma_accuracy*sigmax/dx) - 0.5;
        i2 = ii + fabs(sigma_accuracy*sigmax/dx) + 0.5;
        j1 = jj - fabs(sigma_accuracy*sigmay/dy) - 0.5;
        j2 = jj + fabs(sigma_accuracy*sigmay/dy) + 0.5;
        k1 = kk - fabs(sigma_accuracy*sigmaz/dz) - 0.5;
        k2 = kk + fabs(sigma_accuracy*sigmaz/dz) + 0.5;
        if(i1<0) i1 = 0; if(i2>=X_SIZE) i2 = X_SIZE-1;
        if(j1<0) j1 = 0; if(j2>=Y_SIZE) j2 = Y_SIZE-1;
        if(k1<0) k1 = 0; if(k2>=Z_SIZE) k2 = Z_SIZE-1;
        x0 = x[ii]; y0 = y[jj]; z0 = z[kk];
        value = 0;
        weight = 0;
        for(int i=i1; i<=i2; i++) for(int j=j1; j<=j2; j++) for(int k=k1; k<=k2; k++)
        {
            weight_here = gauss3D(x[i]-x0,y[j]-y0,z[k]-z0, sigmax, sigmay, sigmaz);
            weight += weight_here;
            value  += weight_here*cube(i,j,k);
        }
        cubetemp(ii,jj,kk) = value/weight;
        counter++;
    }
    for(int kk=0;kk<Z_SIZE;kk++) for(int jj=0;jj<Y_SIZE;jj++) for(int ii=0;ii<X_SIZE;ii++) cube(ii,jj,kk) = cubetemp(ii,jj,kk);
}

