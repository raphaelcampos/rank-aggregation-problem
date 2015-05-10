
#include <iostream>
#include <string.h>

using namespace std;
int maior = 0;
int compra = 0;
int vende = 0;

#define max(a,b) (((a)>(b))?(a):(b))

int bf_max_sum_subarray(int A[], int ini, int fim, int &m){

	if(fim - ini + 1 <= 1){
		cout << ini << " - " << fim << endl;
		return A[ini];
	}else{
		//int left = bf_max_sum_subarray(A, ini, fim - 1);
		//int right = bf_max_sum_subarray(A, ini + 1, fim);
		int m = 0;
		return m;
	}
	
}
int sf, si = 0;
int dc_max_sum_subarray(int A[], int ini, int fim){
	if(fim - ini + 1 <= 1){
		return A[ini];
	}else{
		int meio = (ini + fim)/2;
		int max1 = dc_max_sum_subarray(A, ini, meio);
		int max2 = dc_max_sum_subarray(A, meio+1, fim);
		
		int sum;

		int max_1 = sum = A[meio];
		for (int i=meio-1; i >= ini ; i--)
		{
			sum += A[i];
			//max_1 = max(max_1, sum);
			if(max_1<sum){
				max_1 = sum;
				si = i;
			}
		}

		int max_2 = sum = A[meio+1];
		for (int i = meio+2; i <= fim ; i++)
		{
			sum += A[i];
			//max_2 = max(max_2, sum);
			if(max_2<sum){
				max_2 = sum;
				sf = i;
			}
		}

		int max3 = max_2 + max_1;

		int mmm = max(max1,max2);
		return max(mmm,max3);

	}
}



int main(int argc, char const *argv[])
{
	int A[] = {10, 5, -15, 20, 50, -1, 3, -30, 10};

	int m = 0;

	cout << dc_max_sum_subarray(A, 0, 8) << "-" << m << endl;

	cout << si << "-" << sf << endl;


	return 0;
}