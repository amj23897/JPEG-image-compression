#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fcntl.h>
#include <malloc.h>
#include <unistd.h>
#define pi 3.142857

void ComputeDct(int** in,double** out,int x,int y)
{
  int i, j, u, v;
  double s;

  for (i = x; i < x+8; i++)
    for (j = y; j < y+8; j++)
    {
      s = 0;

      for (u = x; u < x+8; u++)
        for (v = y; v < y+8; v++)
          s += in[u][v] * cos((2 * (u-x) + 1) * (i-x) * M_PI / 16) *
                          cos((2 * (v-y) + 1) * (j-y) * M_PI / 16) *
               ((i == x) ? 1 / sqrt(8) : sqrt(2)/sqrt(8)) *
               ((j == y) ? 1 / sqrt(8) : sqrt(2)/sqrt(8));

      out[i][j] = s; // 4;
    }
}

void Compute8x8Idct(double** in, double** out)
{
  int i, j, u, v;
  double s;

  for (i = 0; i < 8; i++)
    for (j = 0; j < 8; j++)
    {
      s = 0;

      for (u = 0; u < 8; u++)
        for (v = 0; v < 8; v++)
          s += in[u][v] * cos((2 * i + 1) * u * M_PI / 16) *
                          cos((2 * j + 1) * v * M_PI / 16) *
               ((u == 0) ? 1 / sqrt(2) : 1.) *
               ((v == 0) ? 1 / sqrt(2) : 1.);

      out[i][j] = s / 4;
    }
}

int* zigZagMatrix(int** arr, int n, int m)
{
    int row = 0, col = 0;
 
    // boolean variable that will 1 if we
    // need to increment 'row' value otherwise
    // 0- if increment 'col' value
    int bool_row_inc = 0;

    int* zig;
    zig=(int *)malloc(m*n*sizeof(int));
 	int cnt=0;
    // Print matrix of lower half zig-zag pattern
    int mn = m<n?m:n;
    for (int len = 1; len <= mn; ++len) {
        for (int i = 0; i < len; ++i) {
            
        	zig[cnt]=arr[row][col];
 			cnt++;

            if (i + 1 == len)
                break;
            // If bool_row_increment value is 1
            // increment row and decrement col
            // else decrement row and increment
            // col
            if (bool_row_inc)
                ++row, --col;
            else
                --row, ++col;
        }
 
        if (len == mn)
            break;
 
        // Update row or col vlaue according
        // to the last increment
        if (bool_row_inc)
            ++row, bool_row_inc = 0;
        else
            ++col, bool_row_inc = 1;
    }
 
    // Update the indexes of row and col variable
    if (row == 0) {
        if (col == m - 1)
            ++row;
        else
            ++col;
        bool_row_inc = 1;
    }
    else {
        if (row == n - 1)
            ++col;
        else
            ++row;
        bool_row_inc = 0;
    }
 
    // Print the next half zig-zag pattern
    int MAX = m>n?m:n - 1;
    for (int len, diag = MAX; diag > 0; --diag) {
 
        if (diag > mn)
            len = mn;
        else
            len = diag;
 
        for (int i = 0; i < len; ++i) {
            zig[cnt]=arr[row][col];
 			cnt++;

            if (i + 1 == len)
                break;
 
            // Update row or col vlaue according
            // to the last increment
            if (bool_row_inc)
                ++row, --col;
            else
                ++col, --row;
        }
 
        // Update the indexes of row and col variable
        if (row == 0 || col == m - 1) {
            if (col == m - 1)
                ++row;
            else
                ++col;
 
            bool_row_inc = 1;
        }
 
        else if (col == 0 || row == n - 1) {
            if (row == n - 1)
                ++col;
            else
                ++row;
 
            bool_row_inc = 0;
        }
    }

    return zig;
}

int* zigZagOrder(int n, int m)
{
    int row = 0, col = 0;
 
    // boolean variable that will 1 if we
    // need to increment 'row' value otherwise
    // 0- if increment 'col' value
    int bool_row_inc = 0;

    int* zig;
    zig=(int *)malloc(m*n*sizeof(int));
 	int cnt=0;
    // Print matrix of lower half zig-zag pattern
    int mn = m<n?m:n;
    for (int len = 1; len <= mn; ++len) {
        for (int i = 0; i < len; ++i) {
            
        	zig[cnt]=row*8+col;
 			cnt++;


            if (i + 1 == len)
                break;
            // If bool_row_increment value is 1
            // increment row and decrement col
            // else decrement row and increment
            // col
            if (bool_row_inc)
                ++row, --col;
            else
                --row, ++col;
        }
 
        if (len == mn)
            break;
 
        // Update row or col vlaue according
        // to the last increment
        if (bool_row_inc)
            ++row, bool_row_inc = 0;
        else
            ++col, bool_row_inc = 1;
    }
 
    // Update the indexes of row and col variable
    if (row == 0) {
        if (col == m - 1)
            ++row;
        else
            ++col;
        bool_row_inc = 1;
    }
    else {
        if (row == n - 1)
            ++col;
        else
            ++row;
        bool_row_inc = 0;
    }
 
    // Print the next half zig-zag pattern
    int MAX = m>n?m:n - 1;
    for (int len, diag = MAX; diag > 0; --diag) {
 
        if (diag > mn)
            len = mn;
        else
            len = diag;
 
        for (int i = 0; i < len; ++i) {
            zig[cnt]=row*8+col;
 			cnt++;

            if (i + 1 == len)
                break;
 
            // Update row or col vlaue according
            // to the last increment
            if (bool_row_inc)
                ++row, --col;
            else
                ++col, --row;
        }
 
        // Update the indexes of row and col variable
        if (row == 0 || col == m - 1) {
            if (col == m - 1)
                ++row;
            else
                ++col;
 
            bool_row_inc = 1;
        }
 
        else if (col == 0 || row == n - 1) {
            if (row == n - 1)
                ++col;
            else
                ++row;
 
            bool_row_inc = 0;
        }
    }

    return zig;
}


void main()
{
	FILE* fp;
    int fp1,fp2,i,j,k,l;
    double a;
    char input_image[100],output_image[100],basis_image[100];
    unsigned char image[3000],image1[3000];
    int npix,nscan,pix,scan,i1,i2;
    printf("Give input image name \n");
    scanf("%s",input_image);
    printf("Give size of image \n");
    scanf("%d%d",&npix,&nscan);
    printf("Input Image is %s npix nscan %d %d \n",input_image,npix,nscan);
    sprintf(output_image,"%s_out",input_image);
    sprintf(basis_image,"%s_basis",input_image);
    fp1= open(input_image,O_RDONLY);
    fp = fopen(basis_image,"w");
    if (fp<0)
    {
        printf("Error in opening %s image \n",basis_image);
        exit(1);
    }
    if (fp1<0)
    {
        printf("Error in opening %s image \n",input_image);
        exit(1);
    }
    int** img;
    img=(int **)malloc(nscan*sizeof(int*));
    for (int i = 0; i < nscan; ++i)
    	img[i]=(int *)malloc(npix*sizeof(int));
    double** dct;
    dct=(double **)malloc(nscan*sizeof(double*));
    for (int i = 0; i < nscan; ++i)
    	dct[i]=(double *)malloc(npix*sizeof(double));
    
    int **imgQuan = (int **)malloc(nscan * sizeof(int *));
    for (i=0; i<nscan; i++)
         imgQuan[i] = (int *)malloc(npix * sizeof(int));
    
    int q[8][8] =     { { 16, 11, 10, 16, 24, 40, 51, 61 },
						{ 12, 12, 14, 19, 26, 58, 60, 55 },
						{ 14, 13, 16, 24, 40, 57, 69, 56 },
						{ 14, 17, 22, 29, 51, 87, 80, 62},
						{ 18, 22, 37, 56, 68, 109, 103, 77 },
						{ 24, 35, 55, 64, 81, 104, 113, 92},
						{ 49, 64, 78, 87, 103, 121, 120, 101 },
						{ 72, 92, 95, 98, 112, 100, 103, 99 } };

	//==========COMPRESSION==========    			
    for (scan=0;scan<nscan;scan++)
    {
        read(fp1,&image[0],npix*sizeof(unsigned char));
        for (pix=0;pix<npix;pix++)
        {
            a=(double)image[pix];
            img[scan][pix] = a-128;
        }
    }
    int** mat88;	//matrix to store 8x8 block
    mat88=(int **)malloc(8*sizeof(int*));
    for (int i = 0; i < 8; ++i)
    	mat88[i]=(int *)malloc(8*sizeof(int));
	    
    int* zig;
    for(scan=0;scan<nscan;scan+=8)
    {
        for(pix=0;pix<npix;pix+=8)
        {
            ComputeDct(img,dct,scan,pix);

            //dividing by quantization matrix
            for (i = scan; i < scan+8; i++) {
                for (j = pix; j < pix+8; j++) {
                    imgQuan[i][j] = dct[i][j]/q[i-scan][j-pix];
                }
            }

            //current 8x8 block
            for (i = scan; i < scan+8; i++) {
                for (j = pix; j < pix+8; j++) {
                    mat88[i-scan][j-pix] = imgQuan[i][j];
                }
            }

            //fuction return 64x1 vector in zig-zag pattern
            zig=zigZagMatrix(mat88,8,8);

            //remove zeros from this vector

            for (i = 63; i >=0; i--)
            {
            	// printf("%d ",zig[i] );
            	if(zig[i]!=0)
            		break;
            }

            //storing length of vector in a file
            fprintf(fp, "%d ",i+1 );
            //storing the values in that file
            for (int k = 0; k < i+1; k++)
            {
				fprintf(fp, "%d ",zig[k] );            	
            }
            fprintf(fp, "\n");


            //this file forms the basis image of jpeg compression     
        }
    }

    //==========DECOMPRESSION==========

    fclose(fp);
    fp = fopen(basis_image,"r");
    if (fp<0)
    {
        printf("Error in opening %s image \n",basis_image);
        exit(1);
    }
    int len,val;
    int* zigdctBlock;
    
    
    int** imgDecomp;
    imgDecomp=(int **)malloc(nscan*sizeof(int *));
    for (int i = 0; i < nscan; ++i)
    	imgDecomp[i]=(int *)malloc(npix*sizeof(int));


    double** dctBlock;
    dctBlock=(double **)malloc(8*sizeof(double*));
    for (int i = 0; i < 8; ++i)
    	dctBlock[i]=(double *)malloc(8*sizeof(double));

    double** idctBlock;
    idctBlock=(double **)malloc(8*sizeof(double*));
    for (int i = 0; i < 8; ++i)
    	idctBlock[i]=(double *)malloc(8*sizeof(double));

    int* order = zigZagOrder(8,8);
    zigdctBlock=(int *)malloc(64*sizeof(int));

  	unsigned char img_rc[npix][nscan];

    for(scan=0;scan<nscan;scan+=8)
    {
        for(pix=0;pix<npix;pix+=8)
        {
			//scanning length
			fscanf(fp,"%d",&len);
		    for (i = 0; i < len; ++i)
		    {
		    	fscanf(fp,"%d",&val );
		    	zigdctBlock[i]=val;
		    	//printf("%d ",val);
		    }
		    //remaining must be zero
		    for(;i<64;i++)
		    {
		    	zigdctBlock[i]=0;
		    }

	        // for reordering of 8*8 matrix
		    for (int i = 0; i < 64; ++i)
	    	{
	    		dctBlock[order[i]/8][order[i]%8]=zigdctBlock[i];
                dctBlock[order[i]/8][order[i]%8]*=q[order[i]/8][order[i]%8];
	    	}	
	     
	        //idct of 8*8 matrix
	        Compute8x8Idct(dctBlock,idctBlock);
	        
	        //converting decompressed image
	        for(i=scan;i<scan+8;i++)
            {
                for(j=pix;j<pix+8;j++)
                {
                  imgDecomp[i][j]=(int)(idctBlock[i-scan][j-pix]+128);
                }
            }
        }
	}

    fp2=creat(output_image,0667);
    if (fp2<0)
	{
	printf("Error in creating output %s image\n",output_image);
	exit(1);
	}

    //saving and converting decompressed image to unsigned char        
    for (scan=0;scan<nscan;scan++)
	{
		for (pix=0;pix<npix;pix++)
		{
		      image1[pix]=(unsigned char)imgDecomp[scan][pix];
		}
	    write(fp2,&image1[0],npix*sizeof(unsigned char));
	}

close(fp1);
close(fp2);
}