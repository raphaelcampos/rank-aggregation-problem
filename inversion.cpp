
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

int sort_and_count(int A[], int ini, int fim){
	if(fim + 1 - ini <= 1){
		return 0;
	}else{
		int meio = (ini + fim)/2;
		return sort_and_count(A, ini, meio) +
		 sort_and_count(A, meio + 1, fim) + 
		 merge_and_count(A, ini, fim, meio +1);
		 ;
	}

}



int main(int argc, char const *argv[])
{
	int A[] = {1, 2, 5,4,3,2,1};

	cout << sort_and_count(A, 0, 6) << endl;

	for (int i = 0; i < 7; ++i)
	{
		cout << " " << A[i];
	}

	return 0;
}