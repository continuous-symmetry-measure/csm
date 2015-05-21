__author__ = 'YAEL'


# def print_output(m, outAtoms, csm, dir, dMin, outFile, localCSM, args_dir, chMinOrder, perm ):
"""
    :param m: molecule dictionary
    :param outAtoms: two dimensional list
    :param csm: float
    :param dir: list of float
    :param dMin: float
    :param outFile: output file (already open)
    :param localCSM: list of float
    :param args_dir: args dictionary
    :param chMinOrder: int
    :param perm: list of int
    :return:
    """

def print_output(output_dict, args_dict):
    """
    :param output_dict:
    :param args_dict:
    :return:
    """

    """
    void printOutput(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out, double* localCSM)
    {

    ////////////////////
        printNumbers
    ////////////////////

        int i, j;
        printf("%s: %.6lf\n", options.opName.c_str(), fabs(csm));
        fprintf(out, "%s: %.4lf\n", options.opName.c_str(), fabs(csm));
        fprintf(out, "SCALING FACTOR: %7lf\n", dMin);

        fprintf(out, "\n INITIAL STRUCTURE COORDINATES\n%i\n", m->size());

    ////////////////////////////
        print initial molecule
    ////////////////////////////

        for (i = 0; i<m->size(); i++){
            fprintf(out, "%3s%10lf %10lf %10lf\n",
                m->symbol(i), m->pos()[i][0], m->pos()[i][1], m->pos()[i][2]);
        }

        for (i = 0; i < m->size(); i++) {
            fprintf(out, "%d ", i + 1);
            for (j = 0; j < m->valency(i); j++) {
                fprintf(out, "%d ", m->adjacent(i, j) + 1);
            }
            fprintf(out, "\n");
        }


    ////////////////////////////////////////////
        print RESULTING STRUCTURE COORDINATES
    ////////////////////////////////////////////

        fprintf(out, "\n RESULTING STRUCTURE COORDINATES\n%i\n", m->size());

        for (i = 0; i<m->size(); i++){
            fprintf(out, "%3s%10lf %10lf %10lf\n",
                m->symbol(i), outAtoms[i][0], outAtoms[i][1], outAtoms[i][2]);
        }

        for (i = 0; i < m->size(); i++) {
            fprintf(out, "%d ", i + 1);
            for (j = 0; j < m->valency(i); j++) {
                fprintf(out, "%d ", m->adjacent(i, j) + 1);
            }
            fprintf(out, "\n");
        }
    /////////////////////////////////////
        print dir
    /////////////////////////////////////

        fprintf(out, "\n DIRECTIONAL COSINES:\n\n");
        fprintf(out, "%lf %lf %lf\n", dir[0], dir[1], dir[2]);


    ////////////////////////////
        print norm
    //////////////////////////////

        if (options.printNorm) {
            printf("NORMALIZATION FACTOR: %7lf\n", m->norm());
            printf("SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf\n", dMin);
            printf("DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);
            printf("NUMBER OF EQUIVALENCE GROUPS: %d\n", m->groupNum());
        }


    /////////////////////////////
        print localCSM
    ////////////////////////////
        if (options.printLocal) {
            double sum = 0;
            fprintf(out, "\nLocal CSM: \n");
            for (i = 0; i < m->size(); i++) {
                sum += localCSM[i];
                fprintf(out, "%s %7lf\n", m->symbol(i), localCSM[i]);
            }
            fprintf(out, "\nsum: %7lf\n", sum);
        }
    }"""
