
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List statint(NumericMatrix K, NumericVector Rtemp, double R0, double r1, double r2) {
    int n = K.nrow();
    List ret;
    
    NumericMatrix Kx(n,n), Ky(n,n), D(n,n), W1(n,n), W2(n,n), W(n,n);
    for (int i = 0; i < n-4; i++) {
        Kx(i,i+1) = 0;
        Ky(i,i+1) = R0 - 2*Rtemp[i+1];
        for (int j = (i+2); j < n-2; j++) {
            double add=0, temp=0;
            NumericVector temp_col;
            temp_col = K(_,j);
            for (int k = (i+1); k < j+1; k++) {
                add +=temp_col[k];
                temp +=Rtemp[k];
            }
            Kx(i,j) = Kx(i,j-1) + 2*add;
            Ky(i,j) = R0 - 2*temp + Kx(i,j);
            D(i,j) = Kx(i,j) - Ky(i,j);
            W1(i,j) = r1*(n-j+i)*Kx(i,j)/n + (j-i)*Ky(i,j)/n;
            W2(i,j) = r2*(n-j+i)*Kx(i,j)/n + (j-i)*Ky(i,j)/n;
            W(i,j) = (n-j+i)*Kx(i,j)/n + (j-i)*Ky(i,j)/n;
        }
    }
    ret["Kx"] = Kx;
    ret["Ky"] = Ky;
    ret["D"] = D;
    ret["W1"] = W1;
    ret["W2"] = W2;
    ret["W"] = W;
    return ret;
}

// [[Rcpp::export]]
NumericVector skew(NumericMatrix K, NumericVector Rtemp, NumericVector Rtemp2, double R0, double R2) {
    int nrow = K.nrow();
    NumericVector out(5), icol, jcol, temp_prod, temp_sum, temp55, temp66, temp77;
    double ijcol;
    
    double temp3=0, temp4=0, temp5=0, temp6=0, temp7=0;
    for (int i = 0; i < nrow-1; i++) {
        icol = K(_,i);
        for (int j = (i+1); j < nrow; j++) {
            double temptemp=0, temp777=0;
            jcol = K(_,j);
            ijcol = icol[j];
            temp_prod = icol*jcol;
            temp_sum = icol + jcol;
            temp77 = 2*(temp_sum*(Rtemp-temp_sum)+temp_prod);
            for (int k = 0; k < nrow; k++) {
                temptemp += temp_prod[k];
                temp777 += temp77[k];
            }
            temp3 += ijcol*temptemp;
            temp4 += ijcol*ijcol*( R0 - (Rtemp[i]+Rtemp[j]-ijcol)*2 );
            temp55 = temp_prod*( Rtemp-temp_sum );
            for (int k = 0; k < nrow; k++) {
                temp5 += temp55[k];
            }
            temp66 = Rtemp - temp_sum;
            temp6 += ijcol*( temp66[i]*temp66[j] - temptemp );
            temp7 += ijcol*( R2-Rtemp2[i]-Rtemp2[j]-temp777+temp77[i]+temp77[j] );
        }
    }
    out[0]=temp3; out[1]=temp4; out[2]=temp5; out[3]=temp6; out[4]=temp7;
    return out;
}
