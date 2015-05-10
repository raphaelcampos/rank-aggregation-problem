
#include <iostream>
#include <string.h>

using namespace std;

int merge_and_count(int A[], int ini, int fim, int meio){

	int *B = (int*)malloc(sizeof(int)*(fim-ini+1));
	int count = 0;

	for (int i = ini; i < meio; ++i)
	{
		B[i] = A[i];
	}
	for (int j = fim; j >= meio; j--)
	{
		B[j] = A[fim - j + meio];
	}

	int i = ini;
	int j = fim;
	
	for (int k = ini; k<=fim; ++k)
	{
		if(B[i] <= B[j]){
			A[k] = B[i];

			i++;
		}else if(B[i] > B[j]){
			A[k] = B[j];
			count += meio - i;
			j--;
		}

	}

	return count;

}

int peak_search(int A[], int ini, int fim){
	if(ini > fim){
		return -1;
	}else{

		int meio = (ini +fim)/2;

		if(A[meio] < A[meio+1] && A[meio] > A[meio-1]){
			return peak_search(A, meio+1, fim);
		}else if(A[meio] > A[meio+1] && A[meio] < A[meio-1]){
			return peak_search(A,ini, meio-1);
		}else if(A[meio] > A[meio+1] && A[meio] > A[meio-1]){
			return meio;
		}else{
			return meio;
		}

	}

}



int main(int argc, char const *argv[])
{
	int A[] = {1, 2, 3, 5, 7, 5};

	cout << A[peak_search(A, 0, 5)] << endl;

	for (int i = 0; i < 7; ++i)
	{
		cout << " " << A[i];
	}

	return 0;
}