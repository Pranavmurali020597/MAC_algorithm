#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define nx 129
#define ny 129
#define Re 400
double u[nx+1][ny+1],u_new[nx+1][ny+1],u_avg[nx][ny];
double v[nx+1][ny+1],v_new[nx+1][ny+1],v_avg[nx][ny];
double p[nx+1][ny+1],p_new[nx+1][ny+1],p_avg[nx][ny];
double RSHU[nx+1][ny+1],RSHV[nx+1][ny+1];
double dx,dy,dt,D[nx+1][ny+1],sf[nx+1][ny+1];
double solver();
double initialisation();
double outputfiles();
int main()
{
	initialisation();
	solver();
	outputfiles();
	return 0;
}
double solver()
{
	int i,j;
	double error,e,en,ed,beta;
	int step;step=0;
	dx=1.0/(double)(nx-1);
	dy=1.0/(double)(ny-1);
	dt=0.001;beta=dx/dy;

    for(i=1;i<=nx;i++)
	{
		u[i][ny+1]=1.0;
		u[i][ny]=1.0;
		u_new[i][ny+1]=1.0;
		u_new[i][ny]=1.0;
	}
	error=1.0;
	while(error>0.000001)
	{   
		step++;
		for(i=2;i<nx;i++)
		{
			for(j=2;j<=ny;j++)
			{
				
					RSHU[i][j]=u[i][j]-
				           (dt/(4*dx))*((u[i][j]+u[i+1][j])*(u[i][j]+u[i+1][j])-
						   (u[i][j]+u[i-1][j])*(u[i][j]+u[i-1][j]))-
						   (dt/(4*dy))*((u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-
						   (u[i][j]+u[i][j-1])*(v[i][j-1]+v[i+1][j-1]))+
						   (dt/(Re*dx*dx))*(u[i-1][j]-2*u[i][j]+u[i+1][j])+
						   (dt/(Re*dy*dy))*(u[i][j-1]-2*u[i][j]+u[i][j+1]);
				
			}
		}
		for(i=2;i<=nx;i++)
		{
			for(j=2;j<ny;j++)
			{
			
					RSHV[i][j]=v[i][j]-(dt/(4*dx))*((v[i][j]+v[i+1][j])*(u[i][j]+u[i][j+1])-
					           (v[i][j]+v[i-1][j])*(u[i-1][j]+u[i-1][j+1]))-
							   (dt/(4*dy))*((v[i][j]+v[i][j+1])*(v[i][j]+v[i][j+1])-
							   (v[i][j]+v[i][j-1])*(v[i][j]+v[i][j-1]))+
							   (dt/(Re*dx*dx))*(v[i-1][j]-2*v[i][j]+v[i+1][j])+
							   (dt/(Re*dy*dy))*(v[i][j-1]-2*v[i][j]+v[i][j+1]);
				
			}
		}
		e=1.0;
		// pressure Gauss-Seidel 
		while(e>0.0001)
		{
			
			for(i=2;i<=nx;i++)
			{
				for(j=2;j<=ny;j++)
				{
					p_new[i][j]=(p_new[i-1][j]+p[i+1][j]+
					(beta*beta)*(p_new[i][j-1]+p[i][j+1])-
					((dx*dx)/dt)*((RSHU[i][j]-RSHU[i-1][j])/dx+
					(RSHV[i][j]-RSHV[i][j-1])/dy))/(2*(1+(beta*beta)));
				}
			}
			for(j=2;j<ny;j++)
			{
				p_new[1][j]=p_new[2][j];
				p_new[nx+1][j]=p_new[nx][j];
			}
			for(i=1;i<=nx;i++)
			{
				p_new[i][1]=p_new[i][2];
				p_new[i][ny+1]=p_new[i][ny];
			}
			for(i=2;i<=nx;i++)
			{
				for(j=2;j<=ny;j++)
				{
					en+=fabs(p[i][j]-p_new[i][j]);
					ed+=fabs(p_new[i][j]);
				}
			}
			e=en/ed;
			for(i=1;i<=(nx+1);i++)
			{
				for(j=1;j<=(ny+1);j++)
				{
					p[i][j]=p_new[i][j];
				}
			}
		} // end of Gauss-Seidel
		
		for(i=2;i<nx;i++)
		{
			for(j=2;j<=ny;j++)
			{
				
				u_new[i][j]=RSHU[i][j]-(dt/dx)*(p_new[i+1][j]-p_new[i][j]);
			}
		}
		for(j=2;j<=(ny);j++)
		{
			u_new[1][j]=0.0;
			u_new[nx][j]=0.0;
		}
		for(i=1;i<=nx;i++)
		{
			u_new[i][1]=-u_new[i][2];
			u_new[i][ny+1]=2-u_new[i][ny];
		}
		for(i=2;i<=nx;i++)
		{
			for(j=2;j<ny;j++)
			{
				v_new[i][j]=RSHV[i][j]-(dt/dy)*(p_new[i][j+1]-p_new[i][j]);
			}
		}
		for(i=1;i<=(nx+1);i++)
		{
			v_new[i][1]=0.0;
			v_new[i][ny]=0.0;
		}
		for(j=2;j<ny;j++)
		{
			
			v_new[1][j]=-v_new[2][j];
			v_new[nx+1][j]=-v_new[nx][j];
		}
		error=0.0;
		for (i=2; i<=nx; i++)
		{
			for (j=2; j<=ny; j++)
			{
				D[i][j] =  ( u_new[i][j]-u[i][j] );
				error = error + fabs(D[i][j]);
			}
		}
		
		if (step%1000 ==1)
		{
	    printf("Error is %5.8lf for the step %d\n", error, step);
		}
		for(i=1;i<=nx;i++)
		{
			for(j=1;j<=(ny+1);j++)
			{
				u[i][j]=u_new[i][j];
			}
		}
		for(i=1;i<=(nx+1);i++)
		{
			for(j=1;j<=ny;j++)
			{
				v[i][j]=v_new[i][j];
			}
		}
	}
	for(j=1;j<ny;j++)
	{
		sf[1][j+1]=sf[i][j]+u[i][j]*dy;	}

	for(j=2;j<=ny;j++)
	{
		for(i=2;i<nx;i++)
		{
			sf[i][j]=sf[i-1][j]+v_new[i][j]*dx;
		}
	}
	for(i=1;i<=nx;i++)
	{
		for(j=1;j<=ny;j++)
		{
			u_avg[i][j]=0.5*(u_new[i][j+1]+u_new[i][j]);
			v_avg[i][j]=0.5*(v_new[i][j]+v_new[i+1][j]);
			p_avg[i][j]=0.25*(p_new[i][j]+p_new[i+1][j]+p_new[i][j+1]+p_new[i+1][j+1]);
		}
	}
}
double initialisation()
{
	printf("*MAC ALGORITHM*");
	int i,j;
	for(i=1;i<=nx;i++)
	{
		for(j=1;j<=(ny+1);j++)
		{
			u[i][j]=0.0;
			u_new[i][j]=0.0;
			
			RSHU[i][j]=0.0;
		}
	}
	for(i=1;i<=(nx+1);i++)
	{
		for(j=1;j<=ny;j++)
		{
			v[i][j]=0.0;
			v_new[i][j]=0.0;
			RSHV[i][j]=0.0;
		}
	}
	for(i=1;i<=(nx+1);i++)
	{
		for(j=1;j<=(ny+1);j++)
		{
			p[i][j]=0.0;
			p_new[i][j]=0.0;
		}
	}
	for(i=1;i<=(nx);i++)
	{
		for(j=1;j<=(ny);j++)
		{
			sf[i][j]=0.0;
		}
	}
}
double outputfiles()
{
	FILE *out,*fout1,*fout2,*fout3;
	dx=1.0/(double)(nx-1);
	dy=1.0/(double)(ny-1);
	double x_axis,y_axis;
	int i,j;
	out=fopen("mac.dat","w");
	fout1=fopen("macuvel.dat","w");
	fout2=fopen("macvvel.dat","w");
	fout3=fopen("macsf.dat","w");
	fprintf( out, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\nZONE F=POINT\n");
	fprintf( out, "I=%d, J=%d\t\n", nx, ny );
	for(i=1;i<=nx;i++)
	{
		x_axis=0.0;
		y_axis=0.0;
		for(j=1;j<=nx;j++)
		{
			x_axis=i-1;x_axis=x_axis*dx;
			y_axis=j-1;y_axis=y_axis*dy;
			fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\n",x_axis,y_axis,u_avg[i][j],v_avg[i][j],p_avg[i][j]);
		}
	}
	// CENTRAL --U
  fprintf(fout1, "VARIABLES=\"U\",\"Y\"\n");
  fprintf(fout1, "ZONE F=POINT\n");
  fprintf(fout1, "I=%d\n", nx);

  for ( j = 1 ; j <=ny ; j++ )
  {
    
    y_axis =  (j-1)*dy;
	fprintf(fout1,"%5.8lf\t%5.8lf\n", u_avg[(nx-1)/2][j], y_axis );

  }// CENTRAL --v
  fprintf(fout2, "VARIABLES=\"X\",\"V\"\n");
  fprintf(fout2, "ZONE F=POINT\n");
  fprintf(fout2, "I=%d\n", nx);
  for (i=1;i<=ny;i++ )
  {
    x_axis =  (i-1)*dx;
	fprintf(fout2,"%5.8lf\t%5.8lf\n",x_axis, v_avg[i][(ny-1)/2]);
  }
  // Stream Function 
  fprintf(fout3, "VARIABLES=\"X\",\"Y\",\"Stream Function\"\n");
  fprintf(fout3, "ZONE F=POINT\n");
  fprintf(fout3, "I=%d,J=%d\n", nx,ny);
  for(j=2;j<ny;j++)
	{
		for(i=1;i<nx;i++)
		{
			x_axis=i-1;x_axis=x_axis*dx;
			y_axis=j-1;y_axis=y_axis*dy;
			fprintf(fout3,"%5.8lf\t%5.8lf\t%5.8lf\n",x_axis,y_axis,sf[i][j]);
		}
	}
	fclose(out);
	fclose(fout1);
	fclose(fout3);
	fclose(fout3);
	return 0;
}
