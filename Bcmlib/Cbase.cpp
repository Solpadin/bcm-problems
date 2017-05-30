#include "stdafx.h"
#include "cbase.h"


/////////////////////////////////////////////////////////
//...подготовка данных для визуализации в формате Surfer;
template <typename T> 
void CDraft<T>::GetSurferFormat(FILE * SURF, FILE * SURF1, FILE * SURF2, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	size_t res;
	if (nd && nd->N > 0 && nd->N1 > 0 && SURF) {
		short int i0 = (short int)nd->N, 
                j0 = (short int)nd->N1;
		T *	 out_F = (T *)new_struct(surfer_dim(type())*sizeof(T)); 
		double X, Y, Z, 
				 min1F = 0., max1F = 1., 
				 min2F = 0., max2F = 1., 
				 min3F = 0., max3F = 1.;
		int hit;
		if (SURF) {
			res = fwrite("DSBB", sizeof(char)*4,  1, SURF);
			res = fwrite(& i0,   sizeof(short int), 1, SURF); res = fwrite(& j0,         sizeof(short int), 1, SURF);
			res = fwrite(nd->X,  sizeof(double),  1, SURF); res = fwrite(nd->X+nd->N-1,  sizeof(double),  1, SURF);
			res = fwrite(nd->Y,  sizeof(double),  1, SURF); res = fwrite(nd->Y+nd->N1-1, sizeof(double),  1, SURF);
			res = fwrite(& min1F, sizeof(double), 1, SURF); res = fwrite(& max1F,        sizeof(double),  1, SURF);
		}
		if (SURF1) {
			res = fwrite("DSBB", sizeof(char)*4,  1, SURF1);
			res = fwrite(& i0,   sizeof(short int), 1, SURF1); res = fwrite(& j0,         sizeof(short int), 1, SURF1);
			res = fwrite(nd->X,  sizeof(double),  1, SURF1); res = fwrite(nd->X+nd->N-1,  sizeof(double),  1, SURF1);
			res = fwrite(nd->Y,  sizeof(double),  1, SURF1); res = fwrite(nd->Y+nd->N1-1, sizeof(double),  1, SURF1);
			res = fwrite(& min2F, sizeof(double), 1, SURF1); res = fwrite(& max2F,        sizeof(double),  1, SURF1);
		}
		if (SURF2) {
			res = fwrite("DSBB", sizeof(char)*4,  1, SURF2);
			res = fwrite(& i0,   sizeof(short int), 1, SURF2); res = fwrite(& j0,         sizeof(short int), 1, SURF2);
			res = fwrite(nd->X,  sizeof(double),  1, SURF2); res = fwrite(nd->X+nd->N-1,  sizeof(double),  1, SURF2);
			res = fwrite(nd->Y,  sizeof(double),  1, SURF2); res = fwrite(nd->Y+nd->N1-1, sizeof(double),  1, SURF2);
			res = fwrite(& min3F, sizeof(double), 1, SURF2); res = fwrite(& max3F,        sizeof(double),  1, SURF2);
		}
		min1F = min2F = min3F = MAX_HIT; 
		max1F = max2F = max3F = MIN_HIT;
      for (int j = 0; j < nd->N1; j++)
      for (int i = 0; i < nd->N;  i++) {
			if (id_axis == AXIS_X) {
				Y = nd->X[i];
				Z = nd->Y[j];
				X = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_Y) {
            Z = nd->X[i];
            X = nd->Y[j];
            Y = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_Z) {
            X = nd->X[i];
            Y = nd->Y[j];
            Z = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_SPH) {
            Z = nd->N2 == 1 ? nd->Z[0] : 1.;
            X = cos(nd->X[i])*sin(nd->Y[j])*Z*1.;
            Y = sin(nd->X[i])*sin(nd->Y[j])*Z*2.;
            Z = cos(nd->Y[j])*Z*3.;
			}
			else
			if (id_axis == AXIS_TIME) {
            X = nd->X[i];
            Y = 0.;
            Z = 0.;
				set_param(6, nd->Y[j]);
			}

			float ff = NOT_HIT; 
			if (nd->hit) hit = nd->hit[i+j*nd->N];
			if (hit != -1)  {
				if (SURF) {
					GetFuncAllValues(X, Y, Z, out_F, hit, _FMF, id_variant, iparam);
               res = fwrite(& (ff = (float)to_double(out_F[0])), sizeof(float), 1, SURF);
               if (min1F > to_double(out_F[0])) min1F = to_double(out_F[0]);
               if (max1F < to_double(out_F[0])) max1F = to_double(out_F[0]);
				}
				if (SURF1) {
               res = fwrite(& (ff = (float)to_double(out_F[1])), sizeof(float), 1, SURF1);
               if (min2F > to_double(out_F[1])) min2F = to_double(out_F[1]);
               if (max2F < to_double(out_F[1])) max2F = to_double(out_F[1]);
				}
				if (SURF2) {
               res = fwrite(& (ff = (float)to_double(out_F[2])), sizeof(float), 1, SURF2);
               if (min3F > to_double(out_F[2])) min3F = to_double(out_F[2]);
               if (max3F < to_double(out_F[2])) max3F = to_double(out_F[2]);
				}
			}
			else {
				if (SURF ) res = fwrite(& ff, sizeof(float), 1, SURF);
				if (SURF1) res = fwrite(& ff, sizeof(float), 1, SURF1);
				if (SURF2) res = fwrite(& ff, sizeof(float), 1, SURF2);
			}
      }
      delete_struct(out_F);

//////////////////////////////////////////////////////////////
//...перезапись максимального и минимального значения функции;
		if (SURF) {
			res = fseek(SURF, sizeof(char)*4+sizeof(short int)*2+sizeof(double)*4, SEEK_SET);
			res = fwrite(& min1F, sizeof(double), 1, SURF);
			res = fwrite(& max1F, sizeof(double), 1, SURF);
			res = fseek(SURF, 0L, SEEK_END);
		}
		if (SURF1) {
			res = fseek(SURF1, sizeof(char)*4+sizeof(short int)*2+sizeof(double)*4, SEEK_SET);
			res = fwrite(& min2F, sizeof(double), 1, SURF1);
			res = fwrite(& max2F, sizeof(double), 1, SURF1);
			res = fseek(SURF1, 0L, SEEK_END);
		}
		if (SURF2) {
			res = fseek(SURF2, sizeof(char)*4+sizeof(short int)*2+sizeof(double)*4, SEEK_SET);
			res = fwrite(& min3F, sizeof(double), 1, SURF2);
			res = fwrite(& max3F, sizeof(double), 1, SURF2);
			res = fseek(SURF2, 0L, SEEK_END);
		}
	}
}


////////////////////////////////
//...строковые варианты функции;
template <typename T> 
void CDraft<T>::GetSurferFormat(char * SURF_FILE, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	char buff[1000]; 
	FILE * SURF = NULL, * SURF1 = NULL, * SURF2 = NULL;
	::strcpy(buff, SURF_FILE); strcat(buff, ".grd");
	if (surfer_dim(type()) > 0) SURF = fopen(buff, "w+b");
	::strcpy(buff, SURF_FILE); strcat(buff, "_1.grd");
	if (surfer_dim(type()) > 1) SURF1 = fopen(buff, "w+b");
	::strcpy(buff, SURF_FILE); strcat(buff, "_2.grd");
	if (surfer_dim(type()) > 2) SURF2 = fopen(buff, "w+b");
	GetSurferFormat(SURF, SURF1, SURF2, nd, _FMF, id_variant, id_axis, iparam);
	if (SURF ) fclose(SURF);
	if (SURF1) fclose(SURF1);
	if (SURF2) fclose(SURF2);
}
template <typename T> 
void CDraft<T>::GetSurferFormat(const char * SURF_FILE, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	char buff[1000]; 
	FILE * SURF = NULL, * SURF1 = NULL, * SURF2 = NULL;
	::strcpy(buff, SURF_FILE); strcat(buff, ".grd");
	if (surfer_dim(type()) > 0) SURF = fopen(buff, "w+b");
	::strcpy(buff, SURF_FILE); strcat(buff, "_1.grd");
	if (surfer_dim(type()) > 1) SURF1 = fopen(buff, "w+b");
	::strcpy(buff, SURF_FILE); strcat(buff, "_2.grd");
	if (surfer_dim(type()) > 2) SURF2 = fopen(buff, "w+b");
	GetSurferFormat(SURF, SURF1, SURF2, nd, _FMF, id_variant, id_axis, iparam);
	if (SURF ) fclose(SURF);
	if (SURF1) fclose(SURF1);
	if (SURF2) fclose(SURF2);
}

//////////////////////////////////////////////////////////
//...подготовка данных для визуализации в формате таблицы;
template <typename T> 
void CDraft<T>::GetDataFormat(FILE * DATA, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	int dim = surfer_dim(type());
	if (nd && nd->N > 0 && nd->N1 > 0 && DATA) {
		T *	 out_F = (T *)new_struct(surfer_dim(type())*sizeof(T));
		double X, Y, Z;
		int hit;
      for (int j = 0; j < nd->N1; j++)
      for (int i = 0; i < nd->N;  i++) {
			if (id_axis == AXIS_X) {
				Y = nd->X[i];
				Z = nd->Y[j];
				X = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_Y) {
            Z = nd->X[i];
            X = nd->Y[j];
            Y = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_Z) {
            X = nd->X[i];
            Y = nd->Y[j];
            Z = nd->N2 == 1 ? nd->Z[0] : 0.;
			}
			else
			if (id_axis == AXIS_SPH) {
            Z = nd->N2 == 1 ? nd->Z[0] : 1.;
            X = cos(nd->X[i])*sin(nd->Y[j])*Z*1.;
            Y = sin(nd->X[i])*sin(nd->Y[j])*Z*2.;
            Z = cos(nd->Y[j])*Z*3.;
			}
			else
			if (id_axis == AXIS_TIME) {
            X = nd->X[i];
            Y = 0.;
            Z = 0.;
				set_param(6, nd->Y[j]);
			}

			float ff = NOT_HIT, f1, f2, f3; 
			if (nd->hit) hit = nd->hit[i+j*nd->N];
			if (hit != -1)  {
				if (DATA) {
					GetFuncAllValues(X, Y, Z, out_F, hit, _FMF, id_variant, iparam);
               f1 = (float)to_double(out_F[0]);
               f2 = (float)to_double(out_F[1]);
               f3 = (float)to_double(out_F[2]);
					if (dim > 2) fprintf(DATA, "%g   %g   %g   %g   %g   %g\n", X, Y, Z, f1, f2, f3); else
					if (dim > 1) fprintf(DATA, "%g   %g   %g   %g   %g\n", X, Y, Z, f1, f2); else
					if (dim > 0) fprintf(DATA, "%g   %g   %g   %g\n", X, Y, Z, f1);
				}
			}
			else {
				if (DATA) {
					if (dim > 2) fprintf(DATA, "%g   %g   %g   %g   %g   %g\n", X, Y, Z, ff, ff, ff); else
					if (dim > 1) fprintf(DATA, "%g   %g   %g   %g   %g\n", X, Y, Z, ff, ff); else
					if (dim > 0) fprintf(DATA, "%g   %g   %g   %g\n", X, Y, Z, ff);
				}
			}
      }
      delete_struct(out_F);
	}
	if (DATA) fclose(DATA);
}

////////////////////////////////
//...строковые варианты функции;
template <typename T> 
void CDraft<T>::GetDataFormat(char * DATA_FILE, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	char buff[1000]; 
	FILE * DATA = NULL;
	::strcpy(buff, DATA_FILE); strcat(buff, ".dat");
	int dim = surfer_dim(type());
	if (dim > 0)   DATA = fopen(buff, "w");
	GetDataFormat (DATA, nd, _FMF, id_variant, id_axis, iparam);
	if (DATA) fclose(DATA);
}
template <typename T> 
void CDraft<T>::GetDataFormat(const char * DATA_FILE, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
	char buff[1000]; 
	FILE * DATA = NULL;
	::strcpy(buff, DATA_FILE); strcat(buff, ".dat");
	int dim = surfer_dim(type());
	if (dim > 0)   DATA = fopen(buff, "w");
	GetDataFormat (DATA, nd, _FMF, id_variant, id_axis, iparam);
	if (DATA) fclose(DATA);
}

/////////////////////////////////////////////////////////////////
//...подготовка данных в формате CSV для пространственной задачи;
template <typename T> 
void CDraft<T>::GetCsvFormat(FILE * CSV, CGrid * nd, int id_variant, int id_centroid, CGrid * bnd)
{
	float ff[9] = {NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT, NOT_HIT}; 
	if (CSV && nd) {
      T * out_F = (T *)new_struct(9*sizeof(T));
      if (out_F && nd->X && nd->Y && nd->Z && nd->geom) {
			switch (id_variant) {
			case 0: fprintf(CSV, ",\"Ux\",\"Uy\",\"Uz\"\n"); break;
			case 1: fprintf(CSV, ",\"txx\",\"tyy\",\"tzz\",\"txy\",\"txz\",\"tyz\"\n"); break;
			}
			Num_Value _FU0 = SPL_VALUE, _FU1 = FLUX_X_VALUE, _FU2 = FLUX_Y_VALUE, _FU3 = FLUX_Z_VALUE;

//////////////////////////////////////////////////////////////////
//...записываем результаты на диск в узлах визуализационной сетки;
          if (! id_centroid && nd->hit)
          for (int k = 0;  k < nd->N; k++) {
				 if (id_variant == 0) {
               GetFuncAllValues(nd->X[k], nd->Y[k], nd->Z[k], out_F, nd->hit[k], _FU0, id_variant);
					fprintf(CSV, "%d,%0.15lg,%0.15lg,%0.15lg\n", k+1/*(nd->hit ? nd->hit[k] : (k+1))*/, 
						ff[0] = (float)to_double(out_F[0]), ff[1] = (float)to_double(out_F[1]), ff[2] = (float)to_double(out_F[2]));
				 }
				 if (id_variant == 1) {
               GetFuncAllValues(nd->X[k], nd->Y[k], nd->Z[k], out_F,   nd->hit[k], _FU1, id_variant);
               GetFuncAllValues(nd->X[k], nd->Y[k], nd->Z[k], out_F+3, nd->hit[k], _FU2, id_variant);
               GetFuncAllValues(nd->X[k], nd->Y[k], nd->Z[k], out_F+6, nd->hit[k], _FU3, id_variant);
               
					fprintf(CSV, "%d,%0.15lg,%0.15lg,%0.15lg,%0.15lg,%0.15lg,%0.15lg\n", k+1/*(nd->hit ? nd->hit[k] : (k+1))*/, 
						ff[0] = (float)to_double(out_F[0]), ff[4] = (float)to_double(out_F[4]), ff[8] = (float)to_double(out_F[8]),
						ff[1] = (float)to_double(out_F[1]), ff[2] = (float)to_double(out_F[2]), ff[5] = (float)to_double(out_F[5]));
				 }
          }
          else

/////////////////////////////////////////////////////////
//...записываем результаты на диск в цетроидах элементов;
          if (! bnd)
          for (int j = 1, k = 0;  k < N; k++) {
               double X0 = 0., Y0 = 0., Z0 = 0., f;
               int id_set_centroid = (B[k].type & ERR_CODE) == ZOOM_BLOCK, l, m, i, num;
//////////////////////
//...centroid setting;
               if (id_set_centroid && B[k].bar->graph) {
                   for (i = num = 0; num < B[k].bar->graph[0]; num++) 
                   if  (1 == B[k].bar->ce[num]->cells_dim()) 
                   for (l = 0;  l < B[k].bar->ce[num]->graph[1];   l++) {
                        X0 += B[k].bar->ce[m = B[k].bar->ce[num]->graph[l+2]]->mp[1];
								Y0 += B[k].bar->ce[m                                ]->mp[2];
                        Z0 += B[k].bar->ce[m                                ]->mp[3];
                        i  += 1;
                   }
                   if (i) {
                       X0 *= (f = 1./i);
                       Y0 *=  f;
                       Z0 *=  f;
                   }
               }
               else {
                   X0 = B[k].mp[1];
                   Y0 = B[k].mp[2];
                   Z0 = B[k].mp[3];
               }
               GetFuncAllValues(X0, Y0, Z0, out_F, k, _FU0, id_variant);
               fprintf(CSV, "%d,%0.15lg,%0.15lg,%0.15lg\n", -nd->geom[j+2], 
						ff[0] = (float)to_double(out_F[0]), ff[1] = (float)to_double(out_F[1]), ff[2] = (float)to_double(out_F[2]));
               j += nd->geom[++j]+1;
          }
          else

///////////////////////////////////////////////////////////////////////////////
//...записываем результаты на диск в цетроидах элементов более подробной сетки;
          if (bnd->geom && 0 < id_centroid && id_centroid <= bnd->geom[0])
			 for (int k = bnd->geom[4*(id_centroid-1)+3]; k < bnd->geom[4*(id_centroid-1)+4]; k++) {
               double X = bnd->X[k], Y = bnd->Y[k], Z = bnd->Z[k];
               int id_block = (int)bnd->hit[k], 
                   ex_block = (int)bnd->nY[k];

               GetFuncAllValues(X, Y, Z, out_F, id_block, _FU0, id_variant);

               fprintf(CSV, "%d,%0.15lg, %0.15lg,%0.15lg\n", ex_block, 
						ff[0] = (float)to_double(out_F[0]), ff[1] = (float)to_double(out_F[1]), ff[2] = (float)to_double(out_F[2]));
          }
		}

///////////////////////
//...addditional point;
      double X0, Y0, Z0;
		int hit = -1;
		X0 = .5;
		Y0 =  1.;
		Z0 = .5;
      Poly_struc_in3D (hit, X0, Y0, Z0);
		GetFuncAllValues(X0, Y0, Z0, out_F, hit, SPL_VALUE, id_variant);
		//fprintf(CSV, "additional point:%0.15lg,%0.15lg,%0.15lg\n", 
		//	ff[0] = (float)to_double(out_F[0]), ff[1] = (float)to_double(out_F[1]), ff[2] = (float)to_double(out_F[2]));
      delete_struct(out_F);
	}
}

////////////////////////////////
//...строковые варианты функции;
template <typename T> 
void CDraft<T>::GetCsvFormat(char * CSV_FILE, CGrid * nd, int id_variant, int id_centroid, CGrid * bnd)
{
	FILE * CSV = fopen(CSV_FILE, "w");
	GetCsvFormat(CSV, nd, id_variant, id_centroid, bnd);
	if (CSV) fclose(CSV);
}
template <typename T> 
void CDraft<T>::GetCsvFormat(const char * CSV_FILE, CGrid * nd, int id_variant, int id_centroid, CGrid * bnd)
{
	FILE * CSV = fopen(CSV_FILE, "w");
	GetCsvFormat(CSV, nd, id_variant, id_centroid, bnd);
	if (CSV) fclose(CSV);
}
