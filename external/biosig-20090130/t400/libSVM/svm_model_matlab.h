const char *model_to_matlab_structure(mxArray *plhs[], int num_of_feature, struct svm_model *model);
const char *matlab_matrix_to_model(struct svm_model *model, const mxArray *matlab_struct);
