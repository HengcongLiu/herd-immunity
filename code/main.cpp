// title: Herd immunity induced by COVID-19 vaccination programs and suppression of epidemics caused by the SARS-CoV-2 Delta variant in China
// author: Hengcong Liu
// email: hc_liu@fudan.edu.cn
// last updated: 2021-11-05
// used to plot the main figure in the main text, including the Re, epidemic curve, and burden, etc.
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

// defined global variants
vector<int> groups;                    // population
vector<double> contraindication;       // rate, contraindication including preganant wemen
vector<int> dose;                      // daily administrated dose
vector<vector<double>> contact_matrix; // contact matrix
vector<vector<double>> susceptibility; // age-specific susceptibility

int **Vfdose;   // VF=individuals with administration of first dose
int *Vsdose_s;  // VS=individuals with administration of second dose but don't gain protection
int **Vsdose_p; // VP=individuals with administration of second dose and would gain protection 14 days later
int *Vsdose_r;  // VR=individuals with administration of second dose and gain protection immediately

// defined functions
vector<vector<string>> readData(string path);                                                                                                            // read files
double NGM(vector<vector<double>> cluster);                                                                                                              // calculate dominant eigenvalue of comtrix
int *doseAlloc(int day, int *Target, int veflag, double rVE, int strategy, int dayItr, int groups_number, double VE, int interval);                      // allocate vaccines
void SLIR(char *outpath, double VE, int veflag, double rVE, int interval, int day0, int strategy, double im, int seed, int susflag, double R0, int sim); // run SLIR model

int main(int argc, char *argv[])
{
    Rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(Rng, 1); // gsl seed

    // variable from input
    char *outpath = argv[1];      // output path
    double VE = atof(argv[2]);    // vaccine efficacy against infection
    int veflag = atoi(argv[3]);   // veflag=1, fVE=0.425/0.507*VE; veflag=2, fVE=0
    double rVE = atof(argv[4]);   // relative vaccine efficacy for individuals aged 3-17 and 60+ compared to adults aged 18-59 years
    int interval = atoi(argv[5]); // time interval between 2 doses
    int day0 = atoi(argv[6]);     // epdemic starting day
    int strategy = atoi(argv[7]); // vaccine strategy
    double im = atof(argv[8]);    // initial fraction of immune population
    int seed = atoi(argv[9]);     // number of seeds
    int susflag = atoi(argv[10]); // susflag=1, heterogeneous susceptibility; susflag=2, homogeneous susceptibility
    double R0 = atof(argv[11]);   // basic reproduction number with/without NPIs

    // read files
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

        SLIR(outpath, VE, veflag, rVE, interval, day0, strategy, im, seed, susflag, R0, sim);
        contact_matrix_tmp.clear();
        contact_matrix.clear();
    }

    groups.clear();
    contraindication.clear();
    susceptibility.clear();
    dose.clear();

    return 0;
}

void SLIR(char *outpath, double VE, int veflag, double rVE, int interval, int day0, int strategy, double im, int seed, int susflag, double R0, int sim)
{
    // output files
    ofstream ofp1(string(outpath) + "/cumulative_coverage_groups.txt", ios::app); // simID_day_coverage3-17_coverage18+_coverage_all
    ofstream ofp2(string(outpath) + "/vaccinated_removed_groups.txt", ios::app);  // simID_day_number of protected individuals by groups
    ofstream ofp3(string(outpath) + "/v_infs.txt", ios::app);                     // simID_cumulative number of vaccinated infected individuals by groups
    ofstream ofp4(string(outpath) + "/nv_infs.txt", ios::app);                    // simID_cumulative number of unvaccinated infected individuals by groups
    ofstream ofp5(string(outpath) + "/infs_immune.txt", ios::app);                // simID_day_newinfs_cuminfs_cumulative vaccine-immune
    ofstream ofp6(string(outpath) + "/daily_R0.txt", ios::app);                   // simID_day_Re

    // local variables
    int group_number = contact_matrix.size();                          // group number
    int day = 0, dayItr = 14, dayChange = 337;                         // simulation day, time to reach max vaccinal protection after vaccinating second dose, time to change vaccination program
    int popa = 0, pop_group1 = 0, pop_group2 = 0;                      // total population, population of 3-11, population of 12+
    int cumdose = 0, cumdose_group1 = 0, cumdose_group2 = 0;           // cumulative number of all/3-11/12+ individuals vaccinating 2 doses
    int newinfs = 0, cuminfs = 0, cumprotected = 0;                    // number of new infecitons/cumulative infections/cumulative vaccinated-recovered
    double eigen = 0, beta = 0, latent = 1.0 / 4.4, gamma = 1.0 / 2.6; // dominant eigenvalue, transmission/latent/recovery rate
    double *sus = new double[group_number];                            // age-specific susceptibility,
    Vfdose = new int *[group_number];                                  // VF=individuals vaccinated first dose
    Vsdose_s = new int[group_number];                                  // VS=individuals vaccinated second dose, and don't gain protection,
    Vsdose_p = new int *[group_number];                                // VP=individuals vaccinated second dose, and gain protection 14 day later
    Vsdose_r = new int[group_number];                                  // VR=individuals vaccinated second dose and gain protection immediately
    int *S = new int[group_number];                                    // susceptible individuals
    int *C = new int[group_number];                                    // indiviudls with contraindications
    int *L = new int[group_number];                                    // latent individuals
    int *I = new int[group_number];                                    // infected individuals
    int *R = new int[group_number];                                    // recovered individuals
    int *iR = new int[group_number];                                   // individuals with natural immunity
    int *VL = new int[group_number];                                   // vaccinated-latent individuals
    int *VI = new int[group_number];                                   // vaccinated-infected individuals
    int *VR = new int[group_number];                                   // vaccinated-recovered individuals
    int *n_sl = new int[group_number];                                 // number of individuals change from S to L
    int *n_cl = new int[group_number];                                 // number of individuals change from C to L
    int *n_li = new int[group_number];                                 // number of individuals change from L to I
    int *n_ir = new int[group_number];                                 // number of individuals change from I to R
    int **n_vfl = new int *[group_number];                             // number of indiciduals change from VF to VL
    int *n_vsl = new int[group_number];                                // number of individuals change from VS to VL
    int **n_vpl = new int *[group_number];                             // number of individuals change from VP to VL
    int *n_vli = new int[group_number];                                // number of individuals change from VL to VI
    int *n_vir = new int[group_number];                                // number of individuals change from VI to VR
    int *Target = new int[group_number];                               // target vaccination population

    if (susflag == 1) // heterogeneous suceptibility
    {
        for (int i = 0; i < group_number; i++)
        {
            int group_index = gsl_rng_uniform_int(Rng, susceptibility[i].size());
            sus[i] = susceptibility[i][group_index];
        }
    }
    else if (susflag == 2) // homogeneous susceptibility
    {
        for (int i = 0; i < group_number; i++)
        {
            sus[i] = 1.0;
        }
    }

    vector<vector<double>> contact_matrix_tmp; // estimated beta
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

    for (int i = 0; i < group_number; i++) // population of 3-11, 12+ and all
    {
        popa += groups[i];
        if (i >= 1 and i <= 3)
        {
            pop_group1 += groups[i];
        }
        else if (i >= 4)
        {
            pop_group2 += groups[i];
        }
    }

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
        Vfdose[i] = new int[interval];
        Vsdose_s[i] = 0;
        Vsdose_p[i] = new int[dayItr];
        Vsdose_r[i] = 0;
        for (int j = 0; j < interval; j++)
        {
            Vfdose[i][j] = 0;
        }
        for (int j = 0; j < dayItr; j++)
        {
            Vsdose_p[i][j] = 0;
        }
    }

    for (int i = 0; i < group_number; i++) // allocate individuals with natural immunity to groups and number is proportional to group population
    {
        int tmpS = S[i] * im, tmpC = C[i] * im;
        S[i] -= tmpS;
        C[i] -= tmpC;
        iR[i] = tmpS + tmpC;
    }

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

    for (day = 1; day <= day0 + 365; day++)
    {
        newinfs = 0;
        if (day >= day0)
        {
            contact_matrix_tmp.clear(); // estimate daily effective reproduction number
            for (int i = 0; i < group_number; i++)
            {
                vector<double> vector_tmp;
                for (int j = 0; j < group_number; j++)
                {
                    vector_tmp.push_back(contact_matrix[i][j] * sus[i] * (groups[i] - L[i] - I[i] - R[i] - VL[i] - VI[i] - VR[i] - Vsdose_r[i] - iR[i]) / groups[j]);
                }
                contact_matrix_tmp.push_back(vector_tmp);
                vector_tmp.clear();
            }
            eigen = NGM(contact_matrix_tmp);
            double Re = beta * eigen / gamma;
            ofp6 << sim << "," << day << "," << Re << endl;

            for (int i = 0; i < group_number; i++) // need to initialize each simulation step
            {
                n_sl[i] = 0;
                n_cl[i] = 0;
                n_li[i] = 0;
                n_ir[i] = 0;
                n_vfl[i] = new int[interval];
                n_vsl[i] = 0;
                n_vpl[i] = new int[dayItr];
                n_vli[i] = 0;
                n_vir[i] = 0;
                for (int j = 0; j < interval; j++)
                {
                    n_vfl[i][j] = 0;
                }
                for (int j = 0; j < dayItr; j++)
                {
                    n_vpl[i][j] = 0;
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

                n_sl[i] = int(round(gsl_ran_binomial(Rng, foi, S[i]))); // change from S to L
                n_cl[i] = int(round(gsl_ran_binomial(Rng, foi, C[i]))); // change from C to L
                for (int j = 0; j < interval; j++)                      // change from VF to VL
                {
                    n_vfl[i][j] = int(round(gsl_ran_binomial(Rng, foi, Vfdose[i][j])));
                }
                n_li[i] = int(round(gsl_ran_binomial(Rng, latent, L[i])));      // change from L to I
                n_ir[i] = int(round(gsl_ran_binomial(Rng, gamma, I[i])));       // change from I to R
                n_vsl[i] = int(round(gsl_ran_binomial(Rng, foi, Vsdose_s[i]))); // change from VS to VL
                for (int j = 0; j < dayItr; j++)                                // change from VP to VL
                {
                    n_vpl[i][j] = int(round(gsl_ran_binomial(Rng, foi, Vsdose_p[i][j])));
                }
                n_vli[i] = int(round(gsl_ran_binomial(Rng, latent, VL[i])));
                n_vir[i] = int(round(gsl_ran_binomial(Rng, gamma, VI[i])));
                newinfs += n_li[i] + n_vli[i];
            }

            for (int i = 0; i < group_number; i++) // update compartments
            {
                int sum_vpl = 0, sum_vfl = 0;
                S[i] -= n_sl[i];
                C[i] -= n_cl[i];
                for (int j = 0; j < interval; j++)
                {
                    Vfdose[i][j] -= n_vfl[i][j];
                    sum_vfl += n_vfl[i][j];
                }
                L[i] += (n_sl[i] + n_cl[i] + sum_vfl - n_li[i]);
                I[i] += (n_li[i] - n_ir[i]);
                R[i] += n_ir[i];
                Vsdose_s[i] -= n_vsl[i];
                for (int j = 0; j < dayItr; j++)
                {
                    Vsdose_p[i][j] -= n_vpl[i][j];
                    sum_vpl += n_vpl[i][j];
                }
                VL[i] += (n_vsl[i] + sum_vpl - n_vli[i]);
                VI[i] += (n_vli[i] - n_vir[i]);
                VR[i] += n_vir[i];
            }
        }
        cuminfs += newinfs;

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

        int *v2dose = doseAlloc(day, Target, veflag, rVE, strategy, dayItr, group_number, VE, interval);

        cumprotected = 0;
        ofp2 << sim << "," << day;
        for (int i = 0; i < group_number; i++)
        {
            ofp2 << "," << Vsdose_r[i];
            cumprotected += Vsdose_r[i];
        }
        ofp2 << endl;

        ofp5 << sim << "," << day << "," << newinfs << "," << cuminfs << "," << cumprotected << endl;

        for (int i = 0; i < group_number; i++)
        {
            cumdose += v2dose[i];
            if (i >= 1 and i <= 3)
            {
                cumdose_group1 += v2dose[i];
            }
            else if (i >= 4)
            {
                cumdose_group2 += v2dose[i];
            }
        }

        ofp1 << sim << "," << day << "," << (double)cumdose_group1 / (double)pop_group1 << "," << (double)cumdose_group2 / (double)pop_group2 << "," << (double)cumdose / (double)popa << endl;

        for (int i = 0; i < group_number; i++) // update S because of vaccinating first dose
        {
            S[i] -= Vfdose[i][0];
        }
    }
    ofp1.close();
    ofp2.close();
    ofp5.close();
    ofp6.close();

    ofp3 << sim;
    ofp4 << sim;
    for (int i = 0; i < group_number; i++)
    {
        ofp3 << "," << (VI[i] + VR[i]);
        ofp4 << "," << (I[i] + R[i]);
    }
    ofp3 << endl;
    ofp4 << endl;
    ofp3.close();
    ofp4.close();
}

int *doseAlloc(int day, int *Target, int veflag, double rVE, int strategy, int dayItr, int groups_number, double VE, int interval)
{
    static int *v2dose = new int[groups_number]; // return variable: represents number of daily new add vaccinated individuals (vaccinating 2 doses) by age groups
    for (int i = 0; i < groups_number; i++)
    {
        v2dose[i] = 0;
    }

    int vaccine_number = 0;                                                                 // number of available dose
    int sum_sdose = 0;                                                                      // number of doses used to vaccinate second dose
    int rest_dose = 0;                                                                      // number of doses used to vaccinate first dose
    int sum1dose = 0;                                                                       // number of individuals who need first dose
    int sum2dose = 0;                                                                       // number of individuals who need second dose
    int **tmp_Vfdose = new int *[groups_number], **tmp_Vsdose_p = new int *[groups_number]; // tmp variables to update VF and VP
    int *rest2dose = new int[groups_number];                                                // rest of individuals who needs but doesn't receive second dose
    for (int i = 0; i < groups_number; i++)
    {
        rest2dose[i] = 0;
        tmp_Vfdose[i] = new int[interval];
        tmp_Vsdose_p[i] = new int[dayItr];
        for (int j = 0; j < interval; j++)
        {
            tmp_Vfdose[i][j] = 0;
        }
        for (int j = 0; j < dayItr; j++)
        {
            tmp_Vsdose_p[i][j] = 0;
        }
    }

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
            sum2dose += Vfdose[i][interval - 1];
        }

        if (vaccine_number < sum2dose) // vaccines are not enough to vaccinate all individuals who need second dose
        {
            int tmp = 0;
            for (int i = 0; i < groups_number - 1; i++)
            {
                double fraction = (double)Vfdose[i][interval - 1] / (double)sum2dose;
                v2dose[i] = int(vaccine_number * fraction);
                tmp += int(vaccine_number * fraction);
            }
            v2dose[groups_number - 1] = vaccine_number - tmp;
            if (Vfdose[groups_number - 1][interval - 1] < vaccine_number - tmp)
            {
                v2dose[groups_number - 1] = Vfdose[groups_number - 1][interval - 1];
            }
            for (int i = 0; i < groups_number; i++)
            {
                rest2dose[i] = Vfdose[i][interval - 1] - v2dose[i];
            }
            for (int i = 0; i < groups_number; i++)
            {
                int p1 = 0, p2 = 0, p3 = 0;
                double nVE = 0;
                if (i <= 4 or i >= 14)
                {
                    nVE = rVE * VE;
                }
                else
                {
                    nVE = VE;
                }
                if (veflag == 1)
                {

                    p1 = int(round(v2dose[i] * (1 - nVE)));                                 // VS
                    p2 = int(round(v2dose[i] * (1 - (double)0.425 / (double)0.507) * nVE)); // VP
                    p3 = v2dose[i] - p1 - p2;                                               // VR
                }
                else if (veflag == 2)
                {
                    p1 = int(round(v2dose[i] * (1 - VE)));
                    p2 = v2dose[i] - p1;
                    p3 = 0;
                }
                Vsdose_s[i] += p1;
                tmp_Vsdose_p[i][0] = p2;
                Vsdose_r[i] += (p3 + Vsdose_p[i][dayItr - 1]);
            }
        }
        else // enough vaccines to vaccinate all individuals who need second dose
        {
            for (int i = 0; i < groups_number; i++)
            {
                v2dose[i] = Vfdose[i][interval - 1];
                sum_sdose += Vfdose[i][interval - 1];
                int p1 = 0, p2 = 0, p3 = 0;
                double nVE = 0;
                if (i <= 4 or i >= 14)
                {
                    nVE = rVE * VE;
                }
                else
                {
                    nVE = VE;
                }
                if (veflag == 1)
                {
                    p1 = int(round(Vfdose[i][interval - 1] * (1 - nVE)));                                 // VS
                    p2 = int(round(Vfdose[i][interval - 1] * (1 - (double)0.425 / (double)0.507) * nVE)); // VP
                    p3 = Vfdose[i][interval - 1] - p1 - p2;                                               // VR
                }
                else if (veflag == 2)
                {
                    p1 = int(round(Vfdose[i][interval - 1] * (1 - VE)));
                    p2 = Vfdose[i][interval - 1] - p1;
                    p3 = 0;
                }
                Vsdose_s[i] += p1;
                tmp_Vsdose_p[i][0] = p2;
                Vsdose_r[i] += (p3 + Vsdose_p[i][dayItr - 1]);
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
                    tmp_Vfdose[i][0] = int(rest_dose * fraction);
                    tmp += tmp_Vfdose[i][0];
                }
                tmp_Vfdose[groups_number - 1][0] = rest_dose - tmp;
                if (tmp_Vfdose[groups_number - 1][0] > Target[groups_number - 1])
                {
                    tmp_Vfdose[groups_number - 1][0] = Target[groups_number - 1];
                }
            }
            else // vaccines are enough
            {
                for (int i = 0; i < groups_number; i++)
                {
                    tmp_Vfdose[i][0] = Target[i];
                }
            }
        }

        // update vaccination stage
        for (int i = 0; i < groups_number; i++)
        {
            for (int j = 0; j < interval - 1; j++)
            {
                tmp_Vfdose[i][j + 1] = Vfdose[i][j];
            }
        }

        for (int i = 0; i < groups_number; i++)
        {
            for (int j = 0; j < dayItr - 1; j++)
            {
                tmp_Vsdose_p[i][j + 1] = Vsdose_p[i][j];
            }
        }

        for (int i = 0; i < groups_number; i++)
        {
            for (int j = 0; j < interval - 1; j++)
            {
                Vfdose[i][j] = tmp_Vfdose[i][j];
            }
            Vfdose[i][interval - 1] = tmp_Vfdose[i][interval - 1] + rest2dose[i];
        }

        for (int i = 0; i < groups_number; i++)
        {
            for (int j = 0; j < dayItr; j++)
            {
                Vsdose_p[i][j] = tmp_Vsdose_p[i][j];
            }
        }
    }

    return v2dose;
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
