
// Convenience function to create a 2D array of zeros
double** zeros(int Nrows, int Ncols) {

    double** array; 
    array = new double*[rows]
    for (int i=0; i<Nrows; i++) {
        array[i] = new double[cols];
    }

    return array;
}

double** createnodes(int Nn)
{
	int M = 3;
    double** nodes = zeros(3, Nn);

    for (int i=0; i<Nn; i++) {
        nodes[i][0] = 3.0;
        nodes[i][1] = 3.0; 
        nodes[i][2] = 0.0;
    }

    return nodes;
}

double** createsources(int Ns)
{

}