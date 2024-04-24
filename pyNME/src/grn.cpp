#include <cstdio>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <map>
#include <mlpack/core.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <boost/math/special_functions/digamma.hpp>

// using namespace std;
using namespace mlpack;



std::vector<std::vector<double>> gen_exp(int n, int m){
    std::vector<std::vector<double>> mat(n, std::vector<double>(m));
    // random_device rd;
    std::mt19937 gen(114514);
    std::uniform_real_distribution<double> dist(0, 10);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            mat[i][j] = dist(gen);
        }
    }
    return mat;
}

double calculateMean(const std::vector<double>& data) {
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / data.size();
}

double calculateStandardDeviation(const std::vector<double>& data, double mean) {
    double squaredSum = 0.0;
    for (const auto& value : data) {
        squaredSum += std::pow(value - mean, 2);
    }
    double variance = squaredSum / data.size();
    return std::sqrt(variance);
}

void kMaximumCorrelation(std::vector<std::vector<double>>& data, int k, std::vector<std::vector<int>>& sortedIndices){
    int nrow = data.size(), ncol = data[0].size();

    std::vector<double> means;
    std::vector<double> stddevs;
    std::vector<std::vector<double>> centered_data(nrow, std::vector<double>(ncol, 0));
    // double centered_data[nrow][ncol];
    std::vector<std::vector<double>> correlation_matrix;

    // Calculate the means and standard deviations for each variable
    for (const auto& row : data) {
        double mean = calculateMean(row);
        double stddev = calculateStandardDeviation(row, mean);
        means.push_back(mean);
        stddevs.push_back(stddev);
    }

    for (int i = 0; i < data.size(); i++){
        std::vector<double> centered_row;
        for (int j = 0; j < data[0].size(); j++){
            centered_data[i][j] = data[i][j] - means[i];
        }
    }

    for (int i = 0; i < nrow; i++){
        std::vector<double> correlation_row;
        for (int j = 0; j < nrow; j++){
            double correlation = 0.0;
            for (int k = 0; k < ncol; k++){
                correlation += (centered_data[i][k] * centered_data[j][k]) / (ncol * stddevs[i] * stddevs[j]);
            }
            correlation_row.push_back(correlation);
        }
        correlation_matrix.push_back(correlation_row);
    }


    std::vector<std::vector<int>> cal_gene(nrow, std::vector<int>(nrow, 0));
    // int row_cnt = 0;
    // std::vector<std::vector<int>> sortedIndices;
    for (const auto& row : correlation_matrix) {
        std::vector<int> indices(row.size());
        std::iota(indices.begin(), indices.end(), 0); // Initialize indices from 0 to row.size()-1

        // Sort indices based on the absolute values in descending order
        std::partial_sort(indices.begin(), indices.begin() + k, indices.end(), [&](size_t i, size_t j) {
            return std::abs(row[i]) > std::abs(row[j]);
        });
        sortedIndices.push_back(indices);      
    }
}


void getSortedIndices(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<std::vector<int>>>& neighborMatrix, int p = 3){
    std::vector<std::vector<std::pair<double, int>>> indexedMatrix;
    std::vector<std::vector<int>> sortedIndices;
    int currRow = 0;
    for (const auto& row: matrix){
        std::vector<std::pair<double, int>> indexedRow;
        for (int i = 0; i < row.size(); i++){
            indexedRow.push_back(std::make_pair(row[i], i));
        }
        std::sort(indexedRow.begin(), indexedRow.end());
        // indexedMatrix.push_back(indexedRow);
        std::vector<int> rowIndices;
        for (const auto& pair : indexedRow) {
            rowIndices.push_back(pair.second);
        }
        // sortedIndices.push_back(rowIndices);
        for (int i = 0; i < row.size(); i++){
            int exp_val = row[rowIndices[i]];
            int p_cnt = 0, left = i - 1, right = i + 1;
            while(p_cnt < p && left >= 0 && right < row.size()){
                int left_val = row[rowIndices[left]], right_val = row[rowIndices[right]];
                if (abs(left_val - exp_val) < abs(right_val - exp_val)){
                    neighborMatrix[currRow][rowIndices[i]][p_cnt] = rowIndices[left];
                    left -= 1;
                }
                else{
                    neighborMatrix[currRow][rowIndices[i]][p_cnt] = rowIndices[right];
                    right += 1;
                }
                p_cnt += 1;
            }
            if (left < 0){
                while (p_cnt < p){
                    neighborMatrix[currRow][rowIndices[i]][p_cnt] = rowIndices[right];
                    right += 1;
                    p_cnt += 1;
                }
            }
            if (right >= row.size()){
                while (p_cnt < p){
                    neighborMatrix[currRow][rowIndices[i]][p_cnt] = rowIndices[left];
                    left -= 1;
                    p_cnt += 1;
                }
            }

        }
        currRow += 1;
    }

    // for (int i = 0; i < 6; i++){
    //     for (int j = 0; j < 5; j++){
    //         for (int k = 0; k < 3; k++){
    //             printf("%d ", neighborMatrix[i][j][k]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n\n");
    // }
}

double MutualInformation(double* vec1, double* vec2, int n_samples, int n_features, int k){
  arma::mat data1(vec1, n_features, n_samples, false);
  arma::mat data2(vec2, n_features, n_samples, false);
  arma::mat data = std::move(arma::join_cols(data1, data2));
  NeighborSearch<NearestNeighborSort, ChebyshevDistance> nn(data);
  arma::Mat<size_t> neighbors;
  arma::mat distances;
  nn.Search(k, neighbors, distances);
  auto half_epsilon = distances.row(k-1);
  arma::mat cheby_dist1(n_samples, n_samples, arma::fill::zeros);
  for (int i = 0; i < n_samples; i++){
    for (int j = i+1; j < n_samples; j++){
      cheby_dist1(i, j) = cheby_dist1(j, i) = (arma::max(abs(data1.col(i) - data1.col(j))) < half_epsilon(i));
    }
  }
  arma::mat cheby_dist2(n_samples, n_samples, arma::fill::zeros);
  for (int i = 0; i < n_samples; i++){
    for (int j = i+1; j < n_samples; j++){
      cheby_dist2(i, j) = cheby_dist2(j, i) = (arma::max(abs(data2.col(i) - data2.col(j))) < half_epsilon(i));
    }
  }

  arma::vec nYNN = arma::sum(cheby_dist1, 1);
  arma::vec nYxNN = arma::sum(cheby_dist2, 1);
  int idN[n_samples];
  memset(idN, 0, sizeof(idN));
  for (int i = 0; i < n_samples; i++){
    idN[i] = ((abs(nYNN(i)+1) > 0.01) && (abs(nYxNN(i)+1) > 0.01));
  }

  double psi_nYNN = 0, psi_nYxNN = 0;
  int sample = 0;
  for (int i = 0; i < n_samples; i++){
    if (idN[i]){
      psi_nYNN += boost::math::digamma(nYNN(i)+1);
      psi_nYxNN += boost::math::digamma(nYxNN(i)+1);
      sample += 1;
    }
  }
  double c = abs(boost::math::digamma(k) - psi_nYNN/sample - psi_nYxNN/sample + boost::math::digamma(n_samples+1)); 
  return c;
}

void read_csv(std::string exp_file, std::vector<std::vector<double>>& exp, std::vector<std::string>& all_genes, std::map<std::string, int>& tf_genes, std::string& species){
    std::ifstream file(exp_file); // Replace "matrix.txt" with the path to your matrix file
    std::vector<std::string> colNames;

    std::string tf_filename;
    if (species == "hs"){
        tf_filename = "hs_tf.txt";
    }
    else if (species == "mu"){
        tf_filename = "mu_tf.txt";
    }
    std::ifstream tf_file(tf_filename);
    std::set<std::string> all_tfs;
    std::string tfs;
    while(std::getline(tf_file, tfs)){
        all_tfs.insert(tfs);
    }
    tf_file.close();

    std::string line;
    if (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string colName;
        while (iss >> colName) {
            colNames.push_back(colName);
        }
    }

    int rowcount = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string rowName;
        iss >> rowName;
        all_genes.push_back(rowName);
        if (all_tfs.find(rowName) != all_tfs.end()){
            tf_genes[rowName] = rowcount;
        }

        std::vector<double> row;
        double value;
        while (iss >> value) {
            row.push_back(value);
        }
        exp.push_back(row);
        rowcount++;
    }

    file.close();
}

void write_csv(const std::vector<std::vector<double>>& matrix,
                      const std::map<std::string, int>& rowNames,
                      const std::vector<std::string>& colNames,
                      const std::string& filename) {
    std::ofstream file(filename);

    // Write column names
    file << ",";
    for (const auto& colName : colNames) {
        file << colName << ",";
    }
    file << "\n";

    // Write matrix rows
    int i = 0;
    for (const auto& pair : rowNames) {
        file << pair.first << ",";
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            file << matrix[i][j] << ",";
        }
        file << "\n";
        i++;
    }

    file.close();
}


void write_csv(const std::vector<std::vector<double>>& matrix,
                      const std::vector<std::string>& rowNames,
                      const std::vector<std::string>& colNames,
                      const std::string& filename) {
    std::ofstream file(filename);

    // Write column names
    file << ",";
    for (const auto& colName : colNames) {
        file << colName << ",";
    }
    file << "\n";

    // Write matrix rows
    for (size_t i = 0; i < matrix.size(); ++i) {
        file << rowNames[i] << ",";
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            file << matrix[i][j] << ",";
        }
        file << "\n";
    }

    file.close();
}

int main(int argc, char** argv){
    // int n_genes = 100;
    // int n_samples = 200;
    int n_neighbors = 12;
    int knn = 5;
    std::string exp_file_name;
    std::string species;
    int computed_genes = 0;
    std::string output_file;
    bool use_all_genes = 0;

    // std::vector<std::vector<double>> exp = gen_exp(n_genes, n_samples);

    // parse argument
    for (int i = 1; i < argc; i++){
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h"){
            std::cout << "Usage:" << std::endl;
            std::cout << "  ./scgrn --cm <int> --mi <int> --exp_file <str> --species {hs|mu} [--computed_genes <int>] [--use_all_genes] --output_grn <str>" << std::endl;
            std::cout << "Help:" << std::endl;
            std::cout << "  --cm: Number of neighbors used for cross-mapping." << std::endl;
            std::cout << "  --mi: Number of neighbors used to calculate mutual information." << std::endl;
            std::cout << "  --exp_file/-f: The imput gene expression file with index as genes and colnames as samples. Must be .tsv file." << std::endl;
            std::cout << "  --species/-s: Species: hs or mu" << std::endl;
            std::cout << "  --computed_genes: Number of target genes to compute." << std::endl;
            std::cout << "  --use_all_genes: Calculate the interactions between all genes if TF imformation is not specified." << std::endl;
            std::cout << "  --output_grn/-o: The directory of output network. Must be .csv file." << std::endl;
            return 0;
        }
        else if (arg == "--cm"){
            n_neighbors = atoi(argv[++i]);
            // std::cout << n_neighbors << std::endl;
        }
        else if (arg == "--mi"){
            knn = atoi(argv[++i]);
            // std::cout << knn << std::endl;
        }
        else if (arg == "--exp_file" || arg == "-f"){
            exp_file_name = argv[++i];
            // std::cout << exp_file_name << std::endl;
        }
        else if (arg == "-s" || arg == "--species"){
            species = argv[++i];
            // std::cout << speices << std::endl;
        }
        else if (arg == "--computed_genes"){
            computed_genes = atoi(argv[++i]);
        }
        else if (arg == "--use_all_genes"){
            use_all_genes = 1;
        }
        else if (arg == "-o" || arg == "--output_grn"){
            output_file = argv[++i];
        }
        else {
            std::cerr << "Unkown arguments: " << arg << std::endl;
            return 1;
        }
    }

    std::vector<std::vector<double>> exp; // Vector to store the matrix
    std::vector<std::string> all_genes;
    std::map<std::string, int> tf_genes;
    // tf_genes: {gene_name: row_pos}
    read_csv(exp_file_name, exp, all_genes, tf_genes, species); // species

    std::cout << "Read files done!" << std::endl;


    int n_genes = exp.size();
    int n_samples = exp[0].size();

    int n_tfs = tf_genes.size();

    std::vector<std::vector<int>> sortedIndices;
    kMaximumCorrelation(exp, computed_genes+1, sortedIndices);
    std::cout << "Calculate correlation matrix done!" << std::endl;



    std::vector<std::vector<std::vector<int>>> neighborMatrix;
    neighborMatrix.resize(n_genes, std::vector<std::vector<int>>(n_samples, std::vector<int>(n_neighbors)));

    getSortedIndices(exp, neighborMatrix);
    std::cout << "Find neighbors done!" << std::endl;

    // double mimat[n_genes][n_genes];
    // memset(mimat, 0, sizeof(mimat));
    std::vector<std::vector<double>> mimat(n_tfs, std::vector<double>(n_genes, 0));
    std::vector<std::string> tf_order;

    clock_t start = clock();

    int row_count = 0;
    for (const auto& pair : tf_genes){
        int i = pair.second;
        tf_order.push_back(pair.first);
        for (int k = 1; k <= computed_genes; k++){
            int j = sortedIndices[i][k];
            double YxNN[n_samples][n_neighbors];
            double YNN[n_samples][n_neighbors];
            memset(YxNN, 0, sizeof(YxNN));
            memset(YNN, 0, sizeof(YNN));
            for (int x = 0; x < n_samples; x++){
                for (int y = 0; y < n_neighbors; y++){
                    YxNN[x][y] = exp[j][neighborMatrix[i][x][y]];
                    YNN[x][y] = exp[j][neighborMatrix[j][x][y]];
                }
            }
            double mi1 = MutualInformation(&YNN[0][0], &YxNN[0][0], n_samples, n_neighbors, knn);
            double mi2 = MutualInformation(&YNN[0][0], &YNN[0][0], n_samples, n_neighbors, knn);
            double mi = mi1/mi2;
            mimat[row_count][j] = mi;
            // printf("%d Done!\n", k);
        }
        row_count++;
        clock_t now = clock();
        double curr = (double)(now - start) / CLOCKS_PER_SEC;
        printf("%d Done in %.2fs!\n", row_count, curr);
    }
    
    write_csv(mimat, tf_order, all_genes, output_file);

    clock_t end = clock();
    double duration = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Running time: %.2f seconds\n", duration);

    // for (int i = 0; i < 5; i++){
    //     for (int j = 0; j < 5; j++){
    //         printf("%.2f ", exp[i][j]);
    //     }
    //     printf("\n");
    // }
    return 0;
}
