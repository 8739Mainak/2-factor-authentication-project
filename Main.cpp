#include<bits/stdc++.h>
using namespace std;
#define PI 3.1415926535898 //rounded upto 13th decimal places
#define MAX 1000000
typedef unsigned long long ll;

int primes[78490], p_cnt = 0; // p_cnt will be 78490 for 1000000 numbers
int threshold = 0; // predetermined threshold value
ll prime_product = 0, seed = 0;

void sieve_primes()
{
    bool is_not_prime[MAX + 1]; // autoinitialization with all false
    is_not_prime[0] = is_not_prime[1] = true;
    for (int j = 2; j <= (MAX / 2); j++) // sieving the multiples of 2;
        is_not_prime[2 * j] = true;
    for (int i = 3; i * i <= MAX; i += 2)
    {
        if (!is_not_prime[i])
        {
            for (int j = i; j <= (MAX / i); j += 2) // sieving the multiples of i;
                is_not_prime[i * j] = 1;
        }
    } // O(nlglgn)
    primes[p_cnt++] = 2;
    for (long int i = 3; i <= MAX; i += 2)
    {
        if (!is_not_prime[i])
            primes[p_cnt++] = i;
    }
}

void generate_seed()
{
    ll p1 = 0, p2 = 0; // two big non-pythagorean primes
    while (p1 < 100 || p1 % 4 != 3)
        p1 = primes[(rand() % p_cnt)];
    while (p2 < 100 || p2 == p1 || p2 % 4 != 3)
        p2 = primes[(rand() % p_cnt)];
    prime_product = p1 * p2;

    while ((seed % p1) == 0 || (seed % p2) == 0) // seed is coprime to prime_product
        seed = ((rand() + rand()) % prime_product) + 2;
}

void assign_matrix(float **matrix, int m, int n)
{
    int mat[m][n];  // randomly generated binary matrix by the seed
    int mat1[m][n];
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            mat[i][j] = seed % 2;
            seed = (seed * seed) % prime_product;
        }
    }

    // Gram-Schmidt orthonormalisation

    for (int i = 0; i < m; i++)
        mat1[i][0] = mat[i][0];
    float num, denum, coeff;
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < j; k++)
        { 
            num = 0.0;
            denum = 0.0;
            for (int l = 0; l < m; l++)
            {
                num = num + mat[l][j] * mat1[l][k];
                denum = denum + mat1[l][k] * mat1[l][k];
            }
            coeff = (denum) ? (num / denum) : 0;
            for (int l = 0; l < m; l++)
                mat1[l][j] = mat[l][j] - coeff * mat1[l][k];
            for (int l = 0; l < m; l++)
                mat[l][j] = mat1[l][j];
        }
        float b = 0.0;
        for (int i = 0; i < m; i++)
            b = b + ((mat1[i][j]) * (mat1[i][j]));
        for (int i = 0; i < m; i++)
            matrix[i][j] = (b) ? (mat1[i][j] / sqrt(b)) : 0;
    }
}
void capture_biometric(float *biometric, int m)
{
    /* In RGB format, every element is a triplet ranged between (0 to 255)
    and they can be replaced by a a distribution that is centered around 0*/
    for (int i = 0; i < m; i++) // randomly generated biometric input
        biometric[i] = (rand() % 256) - 128;
}
void apply_dct(float *dct_of_biometric, float *biometric, int m)
{
    float ci, add, sum;
    for (int i = 0; i < m; i++)
    {
        sum = 0;
        for (int k = 0; k < m; k++)
        {
            add = biometric[k] * cos((2 * k + 1) * i * PI / (2 * m));
            sum += add;
        }
        ci = (i) ? (sqrt(2) / sqrt(m)) : (1 / sqrt(m));
        dct_of_biometric[i] = ci * sum;
    }
}
void generate_biohash(bool *biohash, float **matrix, float *dct_of_biometric, int m, int n)
{
    for (int j = 0; j < n; j++)
    {
        float sum = 0.0;
        for (int i = 0; i < m; i++)
        {
            sum = sum + (matrix[i][j] * dct_of_biometric[i]);
        }
        biohash[j] = (sum > threshold);
    }
}

int main()
{
    printf("\tProject Bio-Hashing : Two factor authentication\n\n\n");
    srand(time(0)); /* time(0) returns the time since Unix timestamp in seconds,
    which will be used as the seed of pseudorandom number generator srand() */

    int m = 5, n = 5; // examplary m*n matrix
    float **matrix = (float **)malloc(m * sizeof(float *));
    for (int i = 0; i < m; i++)
        matrix[i] = (float *)malloc(n * sizeof(float));
    // captured biometric is taken as 1-D matrix for simplicity instead of 2-D
    float *biometric = (float *)malloc(m * sizeof(float));
    float *dct_of_biometric = (float *)malloc(m * sizeof(float));
    bool *biohash = (bool *)malloc(m * sizeof(bool));

    sieve_primes();
    generate_seed();

    cout<<"Seed: "<<seed<<"\t\t"<<"Product of the two prime: "<<prime_product<<"\n\n";

    assign_matrix(matrix, m, n); // seed generated orthonormalized matrix
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            cout<<matrix[i][j]<<" ";
        cout<<"\n";
    }
    cout<<"\n";

    capture_biometric(biometric, m);
    for (int i = 0; i < m; i++)
        cout<<biometric[i]<<" ";
    cout<<"\n";

    apply_dct(dct_of_biometric, biometric, m);
    for (int i = 0; i < m; i++)
        cout<<dct_of_biometric[i]<<" "; 
    cout<<"\n\n";

    generate_biohash(biohash, matrix, dct_of_biometric, m, n);
    for (int i = 0; i < n; i++)
        cout<<biohash[i]<<" ";
    cout<<"\n\n";

    for (int i = 0; i < m; i++)
        free(matrix[i]);
    free(matrix);

    return 0;
}
