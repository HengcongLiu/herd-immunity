// title: Herd immunity induced by COVID-19 vaccination programs and suppression of epidemics caused by the SARS-CoV-2 Delta variant in China
// author: Hengcong Liu
// email: hc_liu@fudan.edu.cn
// last updated: 2021-11-05
// used to calculate the disease burden by using leaky vaccine model
//-----------------------------------------------------------------------------------------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
using namespace std;
gsl_rng *Rng;

vector<int> groups;                    // population
vector<double> contraindication;       // rate, contraindication including preganant wemen
vector<int> dose;                      // daily administrated dose
vector<vector<double>> contact_matrix; // contact matrix
vector<vector<double>> susceptibility; // age-specific susceptibility

int **vaxStu; // individuals with administration of doses, 0-20=first dose, 21-34=second dose, 35=max protection

vector<vector<string>> readData(string path);       // read files
double NGM(vector<vector<double>> cluster);         // calculate dominant eigenvalue of comtrix
int *doseAlloc(int day, int *Target, int strategy); // allocate vaccines
void SLIR(int sim, int strategy);                   // run SLIR model

int main(int argc, char *argv[])
{
    Rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(Rng, 1); // gsl seed

    vector<vector<string>> groups_tmp = readData("../input/population");                 // shanghai population
    vector<vector<string>> contraindication_tmp = readData("../input/contraindication"); // contraindication based on WHO
    vector<vector<string>> susceptibility_tmp = readData("../input/susceptibility");     // age-specific susceptibility
    vector<vector<string>> dose_tmp = readData("../input/dose");                         // daily administrated dose scaled to Shanghai
    for (int i = 0; i < groups_tmp.size(); i++)
    {
        groups.push_back(atoi(groups_tmp[i][1].c_str()));
        contraindication.push_back(atof(contraindication_tmp[i][1].c_str()));
        vector<double> tmp;
        for (int j = 1; j < susceptibility_tmp[i].size(); j++)
        {
            tmp.push_back(atof(susceptibility_tmp[i][j].c_str()));
        }
        susceptibility.push_back(tmp);
        tmp.clear();
    }

    for (int i = 0; i < dose_tmp.size(); i++)
    {
        dose.push_back(atoi(dose_tmp[i][1].c_str()));
    }

    int strategy = 2;
    for (int sim = 0; sim < 100; sim++)
    {
        int sim_index = gsl_rng_uniform_int(Rng, 200); // random choose a contact matrix
        string path_contact_matrix = "../input/shanghai-contact-matrix/contact-matrix" + to_string(sim_index);
        vector<vector<string>> contact_matrix_tmp = readData(path_contact_matrix); // contact matrix for Shanghai
        for (int i = 0; i < contact_matrix_tmp.size(); i++)
        {
            vector<double> tmp;
            for (int j = 1; j < contact_matrix_tmp[i].size(); j++)
            {
                tmp.push_back(atof(contact_matrix_tmp[i][j].c_str()));
            }
            contact_matrix.push_back(tmp);
            tmp.clear();
        }
        SLIR(sim, strategy);
        contact_matrix_tmp.clear();
        contact_matrix.clear();
    }

    groups.clear();
    contraindication.clear();
    susceptibility.clear();
    dose.clear();

    return 0;
}

void SLIR(int sim, int strategy)
{
    ofstream ofp1("../out/v_infs.txt", ios::app);  // simID_cumulative number of vaccinated infected individuals by groups
    ofstream ofp2("../out/nv_infs.txt", ios::app); // simID_cumulative number of unvaccinated infected individuals by groups

    int group_number = contact_matrix.size(), interval = 36, day0 = 367, dayChange = 337;
    double *sus = new double[group_number]; // age-specific susceptibility,
    for (int i = 0; i < group_number; i++)
    {
        int group_index = gsl_rng_uniform_int(Rng, susceptibility[i].size());
        sus[i] = susceptibility[i][group_index];
    }

    int *S = new int[group_number];  // susceptible individuals
    int *C = new int[group_number];  // indiviudls with contraindications
    int *L = new int[group_number];  // latent individuals
    int *I = new int[group_number];  // infected individuals
    int *R = new int[group_number];  // recovered individuals
    int *VL = new int[group_number]; // vaccinated-latent individuals
    int *VI = new int[group_number]; // vaccinated-infected individuals
    int *VR = new int[group_number]; // vaccinated-recovered individuals

    int *n_sl = new int[group_number];     // number of individuals change from S to L
    int *n_cl = new int[group_number];     // number of individuals change from C to L
    int *n_li = new int[group_number];     // number of individuals change from L to I
    int *n_ir = new int[group_number];     // number of individuals change from I to R
    int **n_vsl = new int *[group_number]; // number of indiciduals change from V to VL
    int *n_vli = new int[group_number];    // number of individuals change from VL to VI
    int *n_vir = new int[group_number];    // number of individuals change from VI to VR

    int *Target = new int[group_number]; // target vaccination population
    vaxStu = new int *[group_number];

    double R0 = 6.0, eigen = 0, beta = 0, latent = 1.0 / (6.4 - 2.0), gamma = 1.0 / 2.6, VE = 0.543; // dominant eigenvalue, transmission/latent/recovery rate
    vector<vector<double>> contact_matrix_tmp;                                                       // estimated beta
    for (int i = 0; i < group_number; i++)
    {
        vector<double> vector_tmp;
        for (int j = 0; j < group_number; j++)
        {
            vector_tmp.push_back(contact_matrix[i][j] * sus[i] * groups[i] / groups[j]);
        }
        contact_matrix_tmp.push_back(vector_tmp);
        vector_tmp.clear();
    }
    eigen = NGM(contact_matrix_tmp);
    beta = R0 * gamma / eigen;

    for (int i = 0; i < group_number; i++) // initialize
    {
        S[i] = groups[i] * (1 - contraindication[i]);
        C[i] = groups[i] - S[i];
        L[i] = 0;
        I[i] = 0;
        R[i] = 0;
        VL[i] = 0;
        VI[i] = 0;
        VR[i] = 0;
        vaxStu[i] = new int[interval];
        for (int j = 0; j < interval; j++)
        {
            vaxStu[i][j] = 0;
        }
    }

    int seed = 40;
    for (int i = 0; i < seed; i++) // allocate initial infected individuals
    {
        double random_index = gsl_rng_uniform(Rng);
        if (random_index <= 0.5)
        {
            int group_index = gsl_rng_uniform_int(Rng, group_number);
            I[group_index] += 1;
            S[group_index] -= 1;
        }
        else
        {
            int group_index = gsl_rng_uniform_int(Rng, group_number);
            I[group_index] += 1;
            C[group_index] -= 1;
        }
    }

    for (int day = 1; day <= day0 + 365; day++)
    {
        if (day >= 367)
        {
            for (int i = 0; i < group_number; i++) // need to initialize each simulation step
            {
                n_sl[i] = 0;
                n_cl[i] = 0;
                n_li[i] = 0;
                n_ir[i] = 0;
                n_vsl[i] = new int[36];
                for (int j = 0; j < 36; j++)
                {
                    n_vsl[i][j] = 0;
                }
            }

            for (int i = 0; i < group_number; i++) // estimate number of indiviudals transfer between compartments
            {
                double foi = 0;
                for (int j = 0; j < group_number; j++)
                {
                    double fraction = (double)(I[j] + VI[j]) / (double)groups[j];
                    foi += (beta * sus[i] * contact_matrix[i][j] * fraction);
                }

                n_sl[i] = int(round(gsl_ran_binomial(Rng, foi, S[i])));    // change from S to L
                n_cl[i] = int(round(gsl_ran_binomial(Rng, foi, C[i])));    // change from C to L
                n_li[i] = int(round(gsl_ran_binomial(Rng, latent, L[i]))); // change from L to I
                n_ir[i] = int(round(gsl_ran_binomial(Rng, gamma, I[i])));  // change from I to R
                for (int j = 0; j < 36; j++)
                {
                    if (j <= 20)
                    {
                        n_vsl[i][j] = int(round(gsl_ran_binomial(Rng, foi, vaxStu[i][j])));
                    }
                    else if (j >= 21 and j < 35)
                    {
                        n_vsl[i][j] = int(round(gsl_ran_binomial(Rng, foi * (1 - VE * 0.425 / 0.507), vaxStu[i][j])));
                    }
                    else if (j == 35)
                    {
                        n_vsl[i][j] = int(round(gsl_ran_binomial(Rng, foi * (1 - VE), vaxStu[i][j])));
                    }
                }
                n_vli[i] = int(round(gsl_ran_binomial(Rng, latent, VL[i])));
                n_vir[i] = int(round(gsl_ran_binomial(Rng, gamma, VI[i])));
            }

            for (int i = 0; i < group_number; i++) // update compartments
            {
                int sum_vpl = 0, sum_vfl = 0;
                S[i] -= n_sl[i];
                C[i] -= n_cl[i];
                for (int j = 0; j <= 20; j++)
                {
                    vaxStu[i][j] -= n_vsl[i][j];
                    sum_vfl += n_vsl[i][j];
                }
                L[i] += (n_sl[i] + n_cl[i] + sum_vfl - n_li[i]);
                I[i] += (n_li[i] - n_ir[i]);
                R[i] += n_ir[i];
                for (int j = 21; j < 36; j++)
                {
                    vaxStu[i][j] -= n_vsl[i][j];
                    sum_vpl += n_vsl[i][j];
                }
                VL[i] += (sum_vpl - n_vli[i]);
                VI[i] += (n_vli[i] - n_vir[i]);
                VR[i] += n_vir[i];
            }
        }

        // allocate vaccine
        if (strategy == 1) // target population of 12+ years, then change to 3+ on Nov 1, 2021
        {
            if (day < dayChange)
            {
                for (int i = 0; i < group_number; i++)
                {
                    if (i <= 3)
                    {
                        Target[i] = 0;
                    }
                    else
                    {
                        Target[i] = S[i];
                    }
                }
            }
            else if (day >= dayChange)
            {
                Target[0] = 0;
                for (int i = 1; i < group_number; i++)
                {
                    Target[i] = S[i];
                }
            }
        }
        else if (strategy == 2) // target population of 3+ years
        {
            Target[0] = 0;
            for (int i = 1; i < group_number; i++)
            {
                Target[i] = S[i];
            }
        }
        else if (strategy == 0) // no vaccination program
        {
            for (int i = 0; i < group_number; i++)
            {
                Target[i] = 0;
            }
        }

        int *v1dose = doseAlloc(day, Target, strategy);

        for (int i = 0; i < group_number; i++) // update S because of vaccinating first dose
        {
            S[i] -= v1dose[i];
        }
    }

    ofp1 << sim;
    ofp2 << sim;
    for (int i = 0; i < group_number; i++)
    {
        ofp1 << "," << (VI[i] + VR[i]);
        ofp2 << "," << (I[i] + R[i]);
    }
    ofp1 << endl;
    ofp2 << endl;
    ofp1.close();
    ofp2.close();
}

int *doseAlloc(int day, int *Target, int strategy)
{
    int groups_number = groups.size();
    static int *v1dose = new int[groups_number];
    int *rest2dose = new int[groups_number];
    int *v2dose = new int[groups_number];
    int **tmpStu = new int *[groups_number];
    for (int i = 0; i < groups_number; i++)
    {
        v1dose[i] = 0;
        v2dose[i] = 0;
        rest2dose[i] = 0;
        tmpStu[i] = new int[36];
        for (int j = 0; j < 36; j++)
        {
            tmpStu[i][j] = 0;
        }
    }

    int vaccine_number = 0, sum2dose = 0, sum_sdose = 0, rest_dose = 0, sum1dose = 0;
    if (strategy != 0) // there is vaccination program
    {
        if (day <= dose.size())
        {
            vaccine_number = dose[day - 1];
        }
        else
        {
            vaccine_number = dose[dose.size() - 1];
        }

        for (int i = 0; i < groups_number; i++)
        {
            sum2dose += vaxStu[i][20];
        }

        if (vaccine_number < sum2dose) // vaccines are not enough to vaccinate all individuals who need second dose
        {
            int tmp = 0;
            for (int i = 0; i < groups_number - 1; i++)
            {
                double fraction = (double)vaxStu[i][20] / (double)sum2dose;
                v2dose[i] = int(vaccine_number * fraction);
                tmp += int(vaccine_number * fraction);
            }
            v2dose[groups_number - 1] = vaccine_number - tmp;
            if (vaxStu[groups_number - 1][20] < vaccine_number - tmp)
            {
                v2dose[groups_number - 1] = vaxStu[groups_number - 1][20];
            }
            for (int i = 0; i < groups_number; i++)
            {
                rest2dose[i] = vaxStu[i][20] - v2dose[i];
            }
        }
        else // enough vaccines to vaccinate all individuals who need second dose
        {
            for (int i = 0; i < groups_number; i++)
            {
                v2dose[i] = vaxStu[i][20];
                sum_sdose += vaxStu[i][20];
            }

            rest_dose = vaccine_number - sum_sdose; // vaccinate individuals who needs first dose
            for (int i = 0; i < groups_number; i++)
            {
                sum1dose += Target[i];
            }
            if (sum1dose > rest_dose) // vaccines are not enough to vaccinate all individuals who need first dose
            {
                int tmp = 0;
                for (int i = 0; i < groups_number - 1; i++)
                {
                    double fraction = (double)Target[i] / (double)sum1dose;
                    v1dose[i] = int(rest_dose * fraction);
                    tmp += v1dose[i];
                }
                v1dose[groups_number - 1] = rest_dose - tmp;
                if (rest_dose - tmp > Target[groups_number - 1])
                {
                    v1dose[groups_number - 1] = Target[groups_number - 1];
                }
            }
            else // vaccines are enough
            {
                for (int i = 0; i < groups_number; i++)
                {
                    v1dose[i] = Target[i];
                }
            }
        }

        for (int i = 0; i < groups_number; i++)
        {
            for (int j = 0; j < 35; j++)
            {
                tmpStu[i][j + 1] = vaxStu[i][j];
            }
        }

        for (int i = 0; i < groups_number; i++)
        {
            tmpStu[i][0] = v1dose[i];
            tmpStu[i][21] = v2dose[i];
        }

        for (int i = 0; i < groups_number; i++)
        {
            for (int j = 0; j < 36; j++)
            {
                if (j == 35)
                {
                    vaxStu[i][j] += tmpStu[i][j];
                }
                else if (j == 20)
                {
                    vaxStu[i][j] = tmpStu[i][j] + rest2dose[i];
                }
                else
                {
                    vaxStu[i][j] = tmpStu[i][j];
                }
            }
        }
    }

    return v1dose;
}

vector<vector<string>> readData(string path)
{
    static vector<vector<string>> cluster;
    vector<vector<string>>().swap(cluster);
    vector<string> tmpVector;
    vector<string>().swap(tmpVector);

    char textline[100000];
    char *temp, *next_temp;
    ifstream infile(path);

    while (infile)
    {
        infile.getline(textline, 100000);
        if (infile)
        {
            temp = strtok_r(textline, " ", &next_temp);
            while (temp != NULL)
            {
                tmpVector.push_back(temp);
                temp = strtok_r(NULL, " ", &next_temp);
            }
            cluster.push_back(tmpVector);
            tmpVector.clear();
        }
    }
    infile.close();

    return cluster;
}

double NGM(vector<vector<double>> cluster)
{
    double tmp_return;

    int col = cluster.size();
    double data[col * col];
    int position = 0;

    for (int i = 0; i < col; i++)
    {
        for (int j = 0; j < col; j++)
        {
            data[position] = cluster[i][j];
            position += 1;
        }
    }

    gsl_matrix_view m = gsl_matrix_view_array(data, col, col);
    gsl_vector_complex *eval = gsl_vector_complex_alloc(col);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(col, col);
    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(col);
    gsl_eigen_nonsymmv(&m.matrix, eval, evec, w);
    gsl_eigen_nonsymmv_free(w);
    gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

    gsl_complex eval_0 = gsl_vector_complex_get(eval, 0);
    tmp_return = GSL_REAL(eval_0);

    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);

    return tmp_return;
}
