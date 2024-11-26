#include "Matrix.h"

COOMatrix::COOMatrix(bool is_sym)
{
    this->is_sym = is_sym;
}

COOMatrix::COOMatrix(CSRMatrix&& A)
{
    this->is_sym = A.is_sym;
    this->n_line = A.n_row;
    this->n_col = A.n_col;
    this->nnz = A.nnz;
    this->is_sorted = true;

    for (int i = 0; i < A.n_row; i++) {
        int col_start = A.rowptr[i];
        int col_end = A.rowptr[i + 1];
        for (int j_id = col_start; j_id < col_end; j_id++) {
            this->append(i, A.colind[j_id], A.values[j_id]);
        }
    }
}

void COOMatrix::append(int i, int j, double val)
{
    data.push_back({ i, j, val });
    if (i + 1 > n_line) { n_line = i + 1; }
    if (j + 1 > n_col) { n_col = j + 1; }

    is_sorted = 0;
}

std::tuple<int, int> COOMatrix::shape()
{
    return std::make_tuple(n_line, n_col);
}

std::tuple<int, int, double> COOMatrix::get_data(int index)
{
    return data[index];
}

bool compare_column(const std::tuple<int, int, double>& p, const std::tuple<int, int, double>& q) {
    if (std::get<1>(p) < std::get<1>(q)) return true;
    if (std::get<1>(q) < std::get<1>(p)) return false;
    if (std::get<0>(p) < std::get<0>(q)) return true;
    if (std::get<0>(q) < std::get<0>(p)) return false;
    if (std::get<2>(p) < std::get<2>(q)) return true;
    if (std::get<2>(q) < std::get<2>(p)) return false;
};

void COOMatrix::sum_duplicate(bool col_first)
{
    // Add all duplicate (i, j)
    // Sort current COO Matrix
    sort_tuple(col_first);

    // Temporary data that will be swapped
    std::vector<std::tuple<int, int, double>> temp;

    // Go through current COO Matrix
    int iprev = -1, jprev = -1;
    for (int k = 0; k < data.size(); k++) {
        int icurr = std::get<0>(data[k]);
        int jcurr = std::get<1>(data[k]);
        // If index pair exists, add to previous, else push back new data
        if (icurr != iprev or jcurr != jprev) {
            temp.push_back(data[k]);
        }
        else {
            std::get<2>(temp.back()) += std::get<2>(data[k]);
        }
        iprev = icurr; jprev = jcurr;
    }

    // Swap old data with duplicates to new duplicate free
    data.swap(temp);
    is_sorted = 1;
}

void COOMatrix::sum_duplicate()
{
    // Add all duplicate (i, j)
    // Sort current COO Matrix
    std::sort(data.begin(), data.end());

    // Temporary data that will be swapped
    std::vector<std::tuple<int, int, double>> temp;

    // Go through current COO Matrix
    int iprev = -1, jprev = -1;
    for (int k = 0; k < data.size(); k++) {
        int icurr = std::get<0>(data[k]);
        int jcurr = std::get<1>(data[k]);
        // If index pair exists, add to previous, else push back new data
        if (icurr != iprev or jcurr != jprev) {
            temp.push_back(data[k]);
        }
        else {
            std::get<2>(temp.back()) += std::get<2>(data[k]);
        }
        iprev = icurr; jprev = jcurr;
    }

    // Swap old data with duplicates to new duplicate free
    data.swap(temp);
    is_sorted = 1;
}

void COOMatrix::sort_tuple(bool col_first)
{
    if (col_first) {
        std::sort(data.begin(), data.end(), compare_column);
    }
    else {
        std::sort(data.begin(), data.end());
    }
}

size_t COOMatrix::get_nnz()
{
    if (!is_sorted) {
        // Sort and sum duplicates
        sum_duplicate();
    }
    return data.size();
}

//void COOMatrix::tocsr(CSRMatrix& target_csr)
//{
//    // Make sure matrix is (i,j) sorted and summed duplicates
//    sum_duplicate();
//}

void COOMatrix::print_data()
{
    for (int i = 0; i < data.size(); i++) {
        std::cout << std::get<0>(data[i]) << " " << std::get<1>(data[i]) << " " << std::get<2>(data[i]) << std::endl;
    }
}

void COOMatrix::print_data(int n_line_to_print)
{
    int k = 0;
    int i = std::get<0>(data[0]);
    while (i < n_line_to_print) {
        std::cout << std::get<0>(data[k]) << " " << std::get<1>(data[k]) << " " << std::get<2>(data[k]) << std::endl;
        k++;
        i = std::get<0>(data[k]);
    }
}

void COOMatrix::to_dense(std::vector<std::vector<double>>& Dense)
{
    Dense.resize(n_line, std::vector<double>(n_col));
    sum_duplicate();

    for (int i = 0; i < data.size(); i++) {
        Dense[std::get<0>(data[i])][std::get<1>(data[i])] = std::get<2>(data[i]);
    }
}

// For CSR Matrix

CSRMatrix::CSRMatrix(bool is_sym)
{
    this->is_sym = is_sym;
}

CSRMatrix::CSRMatrix(COOMatrix& m)
{
    // Make sure matrix is (i,j) sorted and summed duplicates
    m.sum_duplicate();

    nnz = m.get_nnz();
    n_row = std::get<0>(m.shape());
    n_col = std::get<1>(m.shape());
    is_sym = m.is_sym;

    colind.resize(nnz);
    rowptr.resize(n_row + 1);
    values.resize(nnz);

    for (int i = 0; i < nnz; i++)
    {
        std::tuple<int, int, double> data = m.get_data(i);
        values[i] = std::get<2>(data);
        colind[i] = std::get<1>(data);
        rowptr[std::get<0>(data) + 1]++;
    }
    for (int i = 0; i < n_row; i++)
    {
        rowptr[i + 1] += rowptr[i];
    }
}

CSRMatrix::CSRMatrix(const CSRMatrix& m)
{
    nnz = m.nnz;
    n_row = m.n_row;
    n_col = m.n_col;
    is_sym = m.is_sym;
    rowptr = m.rowptr;
    colind = m.colind;
    values = m.values;
}

std::tuple<int, int> CSRMatrix::shape()
{
    return std::make_tuple(n_row, n_col);
}

void CSRMatrix::to_dense(std::vector<std::vector<double>>& Dense)
{
    // Make sure target matrix is of correct size
    Dense.resize(n_row, std::vector<double>(n_col));

    // Use row pointers to get column range
    for (int i = 0; i < rowptr.size() - 1; i++) {
        int col_start = rowptr[i];
        int col_end = rowptr[i + 1];

        for (int j = col_start; j < col_end; j++) {
            // Build dense matrix
            Dense[i][colind[j]] = values[j];
        }
    }
}

void CSRMatrix::transpose(CSRMatrix& Res)
{
    Res.n_col = n_col;
    Res.n_row = n_row;
    Res.nnz = nnz;
    Res.is_sym = is_sym;

    CSR_to_CSC(n_col, n_row, rowptr, colind, values, Res.rowptr, Res.colind, Res.values);
}

CSRMatrix& CSRMatrix::operator=(CSRMatrix const& mat)
{
    this->nnz = mat.nnz;
    this->n_row = mat.n_row;
    this->n_col = mat.n_col;
    this->rowptr = mat.rowptr;
    this->colind = mat.colind;
    this->values = mat.values;
    this->is_sym = mat.is_sym;
    return *this;
}

void CSRMatrix::print_data()
{
    for (int i = 0; i < rowptr.size(); i++) {
        std::cout << rowptr[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < colind.size(); i++) {
        std::cout << colind[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < values.size(); i++) {
        std::cout << values[i] << " ";
    }
    std::cout << std::endl;
}

void CSRMatrix::print_data(int n_line)
{
    for (int i = 0; i < n_line; i++) {
        for (int j_id = rowptr[i]; j_id < rowptr[i + 1]; j_id++) {
            std::cout << i << " " << colind[j_id] << " " << values[j_id] << std::endl;
        }
    }
}

void CSRMatrix::store_in_file(std::string filename)
{
    std::ofstream myfile;
    myfile.open(filename + "_val.txt");
    for (int i = 0; i < this->values.size(); i++) {
        myfile << this->values[i] << std::endl;
    }
    myfile.close();
    myfile.open(filename + "_colind.txt");
    for (int i = 0; i < this->colind.size(); i++) {
        myfile << this->colind[i] << std::endl;
    }
    myfile.close();
    myfile.open(filename + "_rowptr.txt");
    for (int i = 0; i < this->rowptr.size(); i++) {
        myfile << this->rowptr[i] << std::endl;
    }
    myfile.close();
}

// For CSC Matrix

CSCMatrix::CSCMatrix(bool is_sym)
{
    this->is_sym = is_sym;
}

CSCMatrix::CSCMatrix(COOMatrix& m)
{
    // Make sure matrix is (i,j) sorted and summed duplicates
    m.sum_duplicate(1);

    nnz = m.get_nnz();
    n_row = std::get<0>(m.shape());
    n_col = std::get<1>(m.shape());
    is_sym = m.is_sym;

    rowind.resize(nnz);
    colptr.resize(n_col + 1);
    values.resize(nnz);

    for (int i = 0; i < nnz; i++)
    {
        std::tuple<int, int, double> data = m.get_data(i);
        values[i] = std::get<2>(data);
        rowind[i] = std::get<0>(data);
        colptr[std::get<1>(data) + 1]++;
    }
    for (int i = 0; i < n_col; i++)
    {
        colptr[i + 1] += colptr[i];
    }
}

CSCMatrix::CSCMatrix(const CSCMatrix& m)
{
    nnz = m.nnz;
    n_row = m.n_row;
    n_col = m.n_col;
    is_sym = m.is_sym;
    rowind = m.rowind;
    colptr = m.colptr;
    values = m.values;
}

std::tuple<int, int> CSCMatrix::shape()
{
    return std::make_tuple(n_row, n_col);
}

void CSCMatrix::to_dense(std::vector<std::vector<double>>& Dense)
{
    // Make sure target matrix is of correct size
    Dense.resize(n_row, std::vector<double>(n_col));

    // Use col pointers to get row range
    for (int j = 0; j < colptr.size() - 1; j++) {
        int row_start = colptr[j];
        int row_end = colptr[j + 1];

        for (int i = row_start; i < row_end; i++) {
            // Build dense matrix
            Dense[rowind[i]][j] = values[i];
        }
    }
}

void CSCMatrix::transpose(CSCMatrix& Res)
{
    Res.n_col = n_col;
    Res.n_row = n_row;
    Res.nnz = nnz;
    Res.is_sym = is_sym;

    CSR_to_CSC(n_col, n_row, colptr, rowind, values, Res.colptr, Res.rowind, Res.values);
}

CSCMatrix& CSCMatrix::operator=(CSCMatrix const& mat)
{
    this->nnz = mat.nnz;
    this->n_row = mat.n_row;
    this->n_col = mat.n_col;
    this->is_sym = mat.is_sym;
    this->rowind = mat.rowind;
    this->colptr = mat.colptr;
    this->values = mat.values;
    return *this;
}

void CSCMatrix::print_data()
{
    for (int i = 0; i < colptr.size(); i++) {
        std::cout << colptr[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < rowind.size(); i++) {
        std::cout << rowind[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < values.size(); i++) {
        std::cout << values[i] << " ";
    }
    std::cout << std::endl;
}

// Operators

template<class binary_op>
void CSR_binop_CSR(CSRMatrix const& A, CSRMatrix const& B, CSRMatrix& C, binary_op op)
{
    // Make sure enough room in result matrix
    int nnzmax = A.nnz + B.nnz;
    int n_row = A.n_row;
    int n_col = A.n_col;

    C.values.resize(nnzmax);
    C.colind.resize(nnzmax);
    C.rowptr.resize(A.rowptr.size());

    int nnz = 0;
    C.rowptr[0] = 0;

    for (int i = 0; i < n_row; i++) {
        int Apos = A.rowptr[i]; int Bpos = B.rowptr[i];
        int Aend = A.rowptr[i + 1]; int Bend = B.rowptr[i + 1];

        while (Apos < Aend && Bpos < Bend) {
            int Aj = A.colind[Apos];
            int Bj = B.colind[Bpos];

            if (Aj == Bj) {
                double result = op(A.values[Apos], B.values[Bpos]);
                if (result != 0) {
                    C.colind[nnz] = Aj;
                    C.values[nnz] = result;
                    nnz++;
                }
                Apos++; Bpos++;
            }
            else if (Aj < Bj) {
                double result = op(A.values[Apos], 0);
                if (result != 0) {
                    C.colind[nnz] = Aj;
                    C.values[nnz] = result;
                    nnz++;
                }
                Apos++;
            }
            else {
                double result = op(0, B.values[Bpos]);
                if (result != 0) {
                    C.colind[nnz] = Bj;
                    C.values[nnz] = result;
                    nnz++;
                }
                Bpos++;
            }
        }
        while (Apos < Aend) {
            double result = op(A.values[Apos], 0);
            if (result != 0) {
                C.colind[nnz] = A.colind[Apos];
                C.values[nnz] = result;
                nnz++;
            }
            Apos++;
        }
        while (Bpos < Bend) {
            double result = op(0, B.values[Bpos]);
            if (result != 0) {
                C.colind[nnz] = B.colind[Bpos];
                C.values[nnz] = result;
                nnz++;
            }
            Bpos++;
        }

        C.rowptr[i + 1] = nnz;
    }

    C.nnz = nnz;
    C.n_row = n_row;
    C.n_col = n_col;

    C.colind.resize(nnz);
    C.values.resize(nnz);
}

CSRMatrix operator-(CSRMatrix& mat1, CSRMatrix& mat2)
{
    auto shape1 = mat1.shape();
    auto shape2 = mat2.shape();
    assert(shape1 == shape2);
    assert(mat1.is_sym == mat2.is_sym);
    CSRMatrix mat_res(mat1.is_sym);

    CSR_binop_CSR(mat1, mat2, mat_res, std::minus<double>());

    return mat_res;
}

CSRMatrix operator-(CSRMatrix& mat1, CSRMatrix&& mat2)
{
    return mat1 - mat2;
}

template<typename T>
CSRMatrix operator*(CSRMatrix& mat, T scal)
{
    CSRMatrix mat_res(mat);

    for (int i = 0; i < mat.nnz; i++)
        mat_res.values[i] *= scal;

    return mat_res;
}

CSRMatrix operator*(double scal, CSRMatrix& mat)
{
    return mat * scal;
}

CSRMatrix operator+(CSRMatrix& mat1, CSRMatrix& mat2)
{
    auto shape1 = mat1.shape();
    auto shape2 = mat2.shape();
    assert(shape1 == shape2);
    assert(mat1.is_sym == mat2.is_sym);
    CSRMatrix mat_res(mat1.is_sym);

    CSR_binop_CSR(mat1, mat2, mat_res, std::plus<double>());

    return mat_res;
}

CSRMatrix operator+(CSRMatrix& mat1, CSRMatrix&& mat2)
{
    return mat1 + mat2;
}

void CSR_to_CSC(int n_row, int n_col, std::vector<int>& Ap, std::vector<int>& Aj, std::vector<double>& Ax\
    , std::vector<int>& Bp, std::vector<int>& Bi, std::vector<double>& Bx)
{
    int nnz = Ap[n_row];

    Bi.resize(nnz, 0);
    Bx.resize(nnz, 0);
    Bp.resize(n_col + 1, 0);

    // Count nnz per column of A
    for (int n = 0; n < nnz; n++) {
        Bp[Aj[n]]++;
    }
    // Cumsum to get colptrs
    for (int col = 0, cumsum = 0; col < n_col; col++) {
        int temp = Bp[col];
        Bp[col] = cumsum;
        cumsum += temp;
    }
    Bp[n_col] = nnz;

    // Filling values of B
    for (int row = 0; row < n_row; row++) {
        for (int jj = Ap[row]; jj < Ap[row + 1]; jj++) {
            int col = Aj[jj];
            int dest = Bp[col];

            Bi[dest] = row;
            Bx[dest] = Ax[jj];

            Bp[col]++;
        }
    }

    // Finalize B
    for (int col = 0, last = 0; col < n_col + 1; col++) {
        int temp = Bp[col];
        Bp[col] = last;
        last = temp;
    }
}

void CSR_to_CSC(CSRMatrix& A, CSCMatrix& B)
{
    B.n_row = A.n_row;
    B.n_col = A.n_col;
    B.nnz = A.nnz;
    B.is_sym = A.is_sym;

    CSR_to_CSC(A.n_row, A.n_col, A.rowptr, A.colind, A.values, B.colptr, B.rowind, B.values);
}

void CSC_to_CSR(CSCMatrix& A, CSRMatrix& B) {
    B.n_row = A.n_row;
    B.n_col = A.n_col;
    B.nnz = A.nnz;
    B.is_sym = A.is_sym;

    CSR_to_CSC(A.n_col, A.n_row, A.colptr, A.rowind, A.values, B.rowptr, B.colind, B.values);
}


//_____________________________________ Vectors ______________________________________

template<class binary_op>
void Vec_binop_Vec(DataVector const& vec1, DataVector const& vec2, DataVector& vec_res, binary_op op)
{
    int length = vec_res.size();
    std::fill(vec_res.begin(), vec_res.end(), 0);

    for (int i = 0; i < length; i++) {
        vec_res[i] = op(vec1[i], vec2[i]);
    }
}

DataVector operator*(CSRMatrix& mat, DataVector& vec)
{
    auto shape_mat = mat.shape();
    int n_row = std::get<0>(shape_mat);
    auto len_vec = vec.size();
    assert(std::get<1>(shape_mat) == len_vec);
    DataVector res_vec(n_row);

    if (mat.is_sym) {
        Vec_mult_Mat_sym(mat, vec, res_vec);
    }
    else {
        Vec_mult_Mat(mat, vec, res_vec);
    }

    return res_vec;
}

DataVector operator*(CSRMatrix&& mat, DataVector& vec)
{
    return mat * vec;
}

void Vec_mult_Mat(CSRMatrix& mat, DataVector& vec, DataVector& res_vec)
{
    for (int i = 0; i < mat.n_row; i++) {
        res_vec[i] = 0;
        for (int j = mat.rowptr[i]; j < mat.rowptr[i + 1]; j++) {
            res_vec[i] += mat.values[j] * vec[mat.colind[j]];
        }
    }
}

void Vec_mult_Mat_sym(CSRMatrix& mat, DataVector& vec, DataVector& res_vec)
{
    for (int i = 0; i < mat.n_row; i++) {
        for (int j_ind = mat.rowptr[i]; j_ind < mat.rowptr[i + 1]; j_ind++) {
            int j = mat.colind[j_ind];
            res_vec[i] += mat.values[j_ind] * vec[j];
            if (i != j) { // Add the symmetric term if not on the diagonal
                res_vec[j] += mat.values[j_ind] * vec[i];
            }
        }
    }
}

DataVector Solve_conditioning(CSRMatrix& L, DataVector& vec)
{
    // Solves system L*L^T*y = vec using forward then backward substitution

    DataVector res1 = Forward_substitution(L, vec);

    CSRMatrix LT;
    L.transpose(LT);
    CSCMatrix LT_CSC;
    CSR_to_CSC(LT, LT_CSC);

    DataVector res2 = Backward_substitution(LT_CSC, vec);

    return res2;
}

DataVector Forward_substitution(CSRMatrix& L, DataVector& vec)
{
    // Solve Matrix system vec = L*y by forward substitution, needing a CSR Matrix

    // Check that inversion makes sense
    assert(vec.size() == L.n_col && L.n_col == L.n_row);
    size_t len = vec.size();
    DataVector res(len);

    for (int i = 0; i < L.n_row; i++)
    {
        double cumsum = 0;
        int row_start = L.rowptr[i];
        int row_end = L.rowptr[i + 1];
        for (int j_ind = row_start; j_ind < row_end - 1; j_ind++) {

            cumsum += L.values[j_ind] * res[L.colind[j_ind]];
        }

        res[i] = (vec[i] - cumsum) / L.values[row_end - 1];
    }

    return res;
}

DataVector Backward_substitution(CSCMatrix& L, DataVector& vec)
{
    // Solve Matrix system vec = L^T*y by backward substitution, needing a CSC Matrix

    // Check that inversion makes sense
    assert(vec.size() == L.n_row && L.n_col == L.n_row);
    size_t len = vec.size();
    DataVector res(len);

    for (int j = L.n_col; j > 0; j--)
    {
        double cumsum = 0;
        int col_start = L.colptr[j - 1];
        int col_end = L.colptr[j];
        for (int i_ind = col_start + 1; i_ind < col_end; i_ind++) {

            cumsum += L.values[i_ind] * res[L.rowind[i_ind]];
        }

        res[j - 1] = (vec[j - 1] - cumsum) / L.values[col_start];
    }

    return res;
}

DataVector operator+(DataVector& vec1, DataVector& vec2)
{
    size_t len1 = vec1.size();
    size_t len2 = vec2.size();
    assert(len1 == len2);

    DataVector res_vec(len1);

    Vec_binop_Vec(vec1, vec2, res_vec, std::plus<double>());

    return res_vec;
}

DataVector operator+(DataVector& vec1, DataVector&& vec2)
{
    return vec1 + vec2;
}

DataVector operator+(DataVector&& vec1, DataVector&& vec2)
{
    return vec1 + vec2;
}

DataVector operator-(DataVector& vec1, DataVector& vec2)
{
    size_t len1 = vec1.size();
    size_t len2 = vec2.size();
    assert(len1 == len2);

    DataVector res_vec(len1);

    Vec_binop_Vec(vec1, vec2, res_vec, std::minus<double>());

    return res_vec;
}

DataVector operator-(DataVector& vec1, DataVector&& vec2)
{
    return vec1 - vec2;
}

DataVector operator*(DataVector& vec, double scal)
{
    DataVector res_vec(vec.size());

    for (int i = 0; i < vec.size(); i++) {
        res_vec[i] = vec[i] * scal;
    }

    return res_vec;
}

DataVector operator*(double scal, DataVector& vec)
{
    return vec * scal;
}

DataVector operator*(double scal, DataVector&& vec)
{
    return vec * scal;
}

double operator^(DataVector& vec1, DataVector& vec2)
{
    auto len1 = vec1.size();
    auto len2 = vec2.size();
    assert(len1 == len2);

    double res = 0;

    for (int i = 0; i < len1; i++) {
        res += vec1[i] * vec2[i];
    }

    return res;
}

double operator^(DataVector& vec1, DataVector&& vec2)
{
    auto len1 = vec1.size();
    auto len2 = vec2.size();
    assert(len1 == len2);

    double res = 0;

    for (int i = 0; i < len1; i++) {
        res += vec1[i] * vec2[i];
    }

    return res;
}



DataVector::DataVector(int n_dof)
{
    data.resize(n_dof);
}

DataVector::DataVector(int n_nodes, int dof_node)
{
    data.resize(n_nodes * dof_node);
    dof_per_node = dof_node;
}

DataVector::DataVector(const DataVector& v)
{
    data = v.data;
    dof_per_node = v.dof_per_node;
}

DataVector::DataVector(std::vector<double>& v)
{
    this->data = v;
    dof_per_node = 1;
}

void DataVector::append(double value)
{
    data.push_back(value);
}

double DataVector::norm()
{
    double sum = 0;
    for (int i = 0; i < data.size(); i++) {
        sum += std::pow(data[i], 2);
    }
    return std::sqrt(sum);
}

void DataVector::print()
{
    for (int i = 0; i < size(); i++) {
        std::cout << data[i] << std::endl;
    }
}

DataVector& DataVector::operator=(DataVector const& v)
{
    this->data = v.data;
    this->dof_per_node = v.dof_per_node;
    return *this;
}

DataVector& DataVector::operator+=(DataVector const& v)
{
    for (int i = 0; i < this->size(); i++) {
        this->data[i] += v[i];
    }
    return *this;
}
