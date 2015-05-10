
#include <iostream>
#include <string.h>

using namespace std;
int maior = 0;
int compra = 0;
int vende = 0;

int sort_and_count(int A[], int ini, int fim, int &compra, int &venda){
	if(fim + 1 - ini <= 2){
		compra = ini;
		venda = fim;

		return A[fim] - A[ini];
	}else{
		int meio = (ini + fim)/2;
		int c1,c2,v1,v2; 
		int m1 = sort_and_count(A, ini, meio, c1,v1);
		int m2 = sort_and_count(A, meio + 1, fim, c2, v2);
	
		if(m1 > m2){ 
			compra = c1;
			venda = v1;
		}else{ 
			m1 = m2;
			compra = c2;
			venda = v2;
		};

		for (int i = ini, j = meio + 1; i <= meio && j <=fim ;)
		{
			if(A[j]-A[i] > m1){
				m1 = A[j]-A[i];
				compra = i;
				venda = j;
			}

			if(A[j]-A[i]>0) j++; else i++;
		}
		return m1;
	}

}



int main(int argc, char const *argv[])
{
	int A[] = {128, 300, 100, 200, 9, 1, 5, 20, 100, 220, 18, 100};

	sort_and_count(A, 0, 9, compra, vende);

	cout << compra << " - " << vende << endl;

	for (int i = 0; i < 7; ++i)
	{
		cout << " " << A[i];
	}

	return 0;
}