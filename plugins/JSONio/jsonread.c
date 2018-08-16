/*
 * $Id: jsonread.c 7190 2017-10-17 11:22:43Z guillaume $
 * Guillaume Flandin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jsmn.h"
#include "mex.h"

/* for MATLAB <= R2008a or GNU Octave */
/*
#if defined(mxSetLogical) || !defined(MATLAB_MEX_FILE)
mxArray *mexCallMATLABWithTrap(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[], const char *functionName) {
    mxArray *mx = NULL;
    const char **fields = (const char *[]){"identifier", "message", "case", "stack"};
    mexSetTrapFlag(1);
    if (mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName)) {
        mx = mxCreateStructMatrix(1, 1, 4, fields);
        mxSetFieldByNumber(mx, 0, 0, mxCreateString("MATLAB:error"));
        mxSetFieldByNumber(mx, 0, 1, mxCreateString(functionName));
        mxSetFieldByNumber(mx, 0, 2, mxCreateCellMatrix(0, 0));
        mxSetFieldByNumber(mx, 0, 3, mxCreateStructMatrix(0, 1, 0, NULL));
        return mx;
    }
    else {
        return NULL;
    }
}
#endif
*/

#if defined(_WIN32) || defined(_WIN64)
  #define strcasecmp _stricmp
#endif

static enum jsonrepsty {
  JSON_REPLACEMENT_STYLE_NOP,
  JSON_REPLACEMENT_STYLE_UNDERSCORE,
  JSON_REPLACEMENT_STYLE_HEX,
  JSON_REPLACEMENT_STYLE_DELETE
} ReplacementStyle;

static int should_convert_to_array(const mxArray *pm) {
    size_t i, j, n, nfields;
    mxClassID cat;
    mwSize ndims, *dim, *d;
    
    if (!mxIsCell(pm)) {
        return 0;
    }
    n = mxGetNumberOfElements(pm);
    if (n) {
        cat = mxGetClassID(mxGetCell(pm, 0));
        ndims = mxGetNumberOfDimensions(mxGetCell(pm, 0));
        dim = (mwSize *)mxGetDimensions(mxGetCell(pm, 0));
    }
    else {
        return 1;
    }
    switch (cat) {
        case mxSTRUCT_CLASS:
        case mxDOUBLE_CLASS:
        case mxLOGICAL_CLASS:
            break;
        default:
            return 0;
    }
    for (i = 1; i < n; i++) {
        if (cat != mxGetClassID(mxGetCell(pm, i))) {
            return 0;
        }
        else if (ndims != mxGetNumberOfDimensions(mxGetCell(pm, i))) {
            return 0;
        }
        d = (mwSize *)mxGetDimensions(mxGetCell(pm, i));
        for (j = 0; j < ndims; j++) {
            if (dim[j] != d[j]) {
                if ( !( (cat == mxDOUBLE_CLASS)
                  && (mxGetNumberOfElements(mxGetCell(pm, i)) <= 1)
                  && (mxGetNumberOfElements(mxGetCell(pm, 0)) <= 1) ) ) {
                    return 0;
                }
            }
        }
        if (cat == mxSTRUCT_CLASS) {
            nfields = mxGetNumberOfFields(mxGetCell(pm, i));
            if (nfields != mxGetNumberOfFields(mxGetCell(pm, 0))) {
                return 0;
            }
            for (j = 0; j < nfields; j++) {
                if (mxGetFieldNumber(mxGetCell(pm, 0), mxGetFieldNameByNumber(mxGetCell(pm, i), j)) != j) {
                    return 0;
                }
            }
        }
    }
    if (cat == mxDOUBLE_CLASS) {
        for (i = 0; i < n; i++) {
            if (mxGetNumberOfElements(mxGetCell(pm, i)) == 0) {
                mxDestroyArray(mxGetCell(pm, i));
                mxSetCell((mxArray *)pm, i, mxCreateDoubleScalar(mxGetNaN()));
            }
        }
    }
    return 1;
}

static int setup_for_cell2mat(mxArray *pm) {
    size_t i, n, a, b;
    mxClassID cat;
    mxArray *cell;
    mwSize d, *olddim = NULL, *newdim = NULL;
    
    n = mxGetNumberOfElements(pm);
    if (!n) {
        return 0;
    }
    cell = mxGetCell(pm, 0);
    cat  = mxGetClassID(cell);
    if ((cat != mxDOUBLE_CLASS) && (cat != mxLOGICAL_CLASS)) {
        return 0;
    }
    if ((mxGetNumberOfDimensions(cell) == 2) && (mxGetM(cell) == 1)) {
        return 0;
    }
    if ((mxGetNumberOfDimensions(cell) == 2) && (mxGetN(cell) == 1)) {
        mxSetN(pm, mxGetM(pm));
        mxSetM(pm, 1);
    }
    else {
        d = mxGetNumberOfDimensions(pm);
        olddim = (mwSize *)mxGetDimensions(pm);
        newdim = mxMalloc((d+1) * sizeof(*newdim));
        for (i=0;i<d;i++) {
            newdim[i] = 1;
        }
        newdim[d] = olddim[0];
        mxSetDimensions(pm, newdim, d+1);
        /*mxFree(newdim);*/
    }
    return 1; /* cell2mat's output will require a call to permute */
}

static int create_struct(char *js, jsmntok_t *tok, mxArray **mx);

static char * get_string(char *js, int start, int end) {
    int i, j;
    for (i = start, j = start; j < end; i++, j++) {
        if (js[j] == '\\') {
            switch (js[j+1]) {
                case '"':
                case '\\':
                case '/':
                    js[i] = js[++j];
                    break;
                case 'b':
                    js[i] = '\b'; j++;
                    break;
                case 'f':
                    js[i] = '\f'; j++;
                    break;
                case 'n':
                    js[i] = '\n'; j++;
                    break;
                case 'r':
                    js[i] = '\r'; j++;
                    break;
                case 't':
                    js[i] = '\t'; j++;
                    break;
                case 'u':
                    /* four-hex-digits */
                    js[i] = 'X'; j+=5;
                    break;
                default:
                    mexPrintf("Unexpected backslash escape.");
            }
        }
        else {
            if (i != j) js[i] = js[j];
        }
    }
    js[i] = '\0';
    js[end] = '\0'; /* not necessary sensu stricto */
    return js + start;
}

static char * valid_fieldname_underscore(char *field) {
    /* matlab.lang.makeValidName
       with ReplacementStyle == 'underscore' and Prefix == 'x' */
    char *re = field, *wr = field;
    int beg = 1;
    while (re[0] != '\0') {
        if ((re[0] == 32) || (re[0] == 9)) { /* ' ' or \t */
            if (re[1] == '\0') {
                break;
            }
            else if ((beg != 1) && (re[1] >= 'a') && (re[1] <= 'z')) {
                re[1] -= 32; /* 'a'-'A' */
            }
            re++;
        }
        else {
            if ( ((re[0] >= '0') && (re[0] <= '9')) 
              || ((re[0] >= 'a') && (re[0] <= 'z'))
              || ((re[0] >= 'A') && (re[0] <= 'Z'))) {
                if (re != wr) { wr[0] = re[0]; }
                wr++;
            }
            else {
                if (ReplacementStyle == JSON_REPLACEMENT_STYLE_UNDERSCORE) {
                    wr[0] = '_'; wr++;
                }
            }
            re++;
            beg = 0;
        }
    }
    wr[0] = '\0';
    if ( (field[0] == '\0') || (field[0] == '_')
      || ((field[0] >= '0') && (field[0] <= '9')) ) {
        field = field - 1;
        field[0] = 'x';
    }
    return field;
}

static char * valid_fieldname_hex(char *field) {
    /* matlab.lang.makeValidName
       with ReplacementStyle == 'hex' and Prefix == 'x' */
    char *re = field, *wr = NULL, *ret = NULL, *str = NULL;
    int sts, beg = 1;
    size_t len = 4*strlen(field)+2; /* 'x' + all '0x??' + \0 */
    mxArray *mx = NULL, *ma = NULL;
    ret = (char *)malloc(len);
    wr = ret;
    if (! (((re[0] >= 'a') && (re[0] <= 'z'))
        || ((re[0] >= 'A') && (re[0] <= 'Z'))) ) {
        wr[0] = 'x'; wr++;
    }
    while (re[0] != '\0') {
        if ((re[0] == 32) || (re[0] == 9)) { /* ' ' or \t */
            if (re[1] == '\0') {
                break;
            }
            else if ((beg != 1) && (re[1] >= 'a') && (re[1] <= 'z')) {
                re[1] -= 32; /* 'a'-'A' */
            }
            re++;
        }
        else {
            if ( ((re[0] >= 'a') && (re[0] <= 'z'))
              || ((re[0] >= 'A') && (re[0] <= 'Z'))
              || ((((re[0] >= '0') && (re[0] <= '9'))
              || (re[0] == '_')) && !beg) ) {
                wr[0] = re[0]; wr++;
            }
            else {
                /* could also use sprintf(wr,"%02X",re[0])
                   but might call unicode2native in the future */
                wr[0] = '0'; wr++;
                wr[0] = 'x'; wr++;
                wr[0] = re[0];
                wr[1] = '\0';
                ma = mxCreateString((const char *)wr);
                sts = mexCallMATLAB(1, &mx, 1, &ma, "dec2hex");
                if (sts != 0) {
                    mexErrMsgTxt("Cannot convert to hexadecimal representation.");
                }
                mxDestroyArray(ma);
                str = mxArrayToString(mx);
                mxDestroyArray(mx);
                wr[0] = str[0]; wr++;
                wr[0] = str[1]; wr++;
                mxFree(str);
            }
            re++;
            beg = 0;
        }
    }
    wr[0] = '\0';
    return ret;
}

static char * valid_fieldname(char *field, int *need_free) {
    switch (ReplacementStyle) {
        case JSON_REPLACEMENT_STYLE_NOP:
            /* nop */
            *need_free = 0;
            break;
        case JSON_REPLACEMENT_STYLE_UNDERSCORE:
        case JSON_REPLACEMENT_STYLE_DELETE:
            field = valid_fieldname_underscore(field);
            *need_free = 0;
            break;
        case JSON_REPLACEMENT_STYLE_HEX:
            field = valid_fieldname_hex(field);
            *need_free = 1;
            break;
        default:
            mexErrMsgTxt("Unknown ReplacementStyle.");
            break;
    }
    return field;
}

static int primitive(char *js, jsmntok_t *tok, mxArray **mx) {
    mxArray *ma = NULL;
    int sts;
    switch (js[tok->start]) {
        case 't' :
            *mx =  mxCreateLogicalScalar(1);
            break;
        case 'f' :
            *mx =  mxCreateLogicalScalar(0);
            break;
        case 'n' :
            *mx =  mxCreateDoubleMatrix(0,0,mxREAL);
            break;
        default: /* '-', '0'..'9' */
            ma =  mxCreateString(get_string(js, tok->start, tok->end));
            sts = mexCallMATLAB(1, mx, 1, &ma, "str2double");
            if (sts != 0) {
                mexErrMsgTxt("Conversion from string to double failed.");
            }
            mxDestroyArray(ma);
            break;
    }
    return 1;
}

static int value(char *js, jsmntok_t *tok, mxArray **mx) {
    *mx = mxCreateString(get_string(js, tok->start, tok->end));
    return 1;
}

static int array(char *js, jsmntok_t *tok, mxArray **mx) {
    int i, j;
    mxArray *ma = NULL;
    mxArray *array[1], *parray[2];
    mwSize d;
    int perm, sts;
    
    *mx = mxCreateCellMatrix(tok->size, 1);
    for (i = 0, j = 0; i < tok->size; i++) {
        j += create_struct(js, tok+1+j, &ma);
        mxSetCell(*mx, i, ma);
    }
    
    /* Convert cell array into array when required */
    if (should_convert_to_array(*mx)) {
        perm = setup_for_cell2mat(*mx);
        sts = mexCallMATLAB(1, array, 1, mx, "cell2mat");
        if (sts == 0) {
            mxDestroyArray(*mx);
            *mx = *array;
            
            if (perm) {
                d = mxGetNumberOfDimensions(*mx);
                parray[0] = *mx;
                parray[1] = mxCreateNumericMatrix(1, d, mxDOUBLE_CLASS, mxREAL);
                mxGetPr(parray[1])[0] = d;
                for (i=1;i<d;i++) {
                    mxGetPr(parray[1])[i] = i;
                }
                sts = mexCallMATLAB(1, mx, 2, parray, "permute");
                mxDestroyArray(parray[1]);
                if (sts != 0) {
                    mexErrMsgTxt("Call to permute() failed.");
                }
            }
        }
        else {
            mexPrintf("Call to cell2mat() failed.\n");
        }
    }
    return j+1;
}

static int object(char *js, jsmntok_t *tok, mxArray **mx) {
    int i, j, k, need_free = 0;
    mxArray *ma = NULL;
    char *field = NULL;
    if (tok->size == 0) {
        *mx = mxCreateStructMatrix(1, 1, 0, NULL);
        return 1;
    }
    for (i = 0, j = 0, k = 0; i < tok->size; i++) {
        field = get_string(js, (tok+1+j)->start, (tok+1+j)->end); /* check it is a JSMN_STRING */
        field = valid_fieldname(field, &need_free);
        j++;
        if (i == 0) {
            *mx = mxCreateStructMatrix(1, 1, 1, (const char**)&field);
        }
        else {
            k = mxGetFieldNumber(*mx, field);
            if (k != -1) {
                mexWarnMsgTxt("Duplicate key.");
                ma = mxGetFieldByNumber(*mx, 0, k);
                mxRemoveField(*mx, k);
                mxDestroyArray(ma);
            }
            k = mxAddField(*mx, field);
            if (k == -1)
                mexErrMsgTxt("mxAddField()");
        }
        if (need_free) { free(field); }
        j += create_struct(js, tok+1+j, &ma);
        mxSetFieldByNumber(*mx, 0, k, ma);
    }
    return j+1;
}

static int create_struct(char *js, jsmntok_t *tok, mxArray **mx) {
    if (tok->type == JSMN_PRIMITIVE) {
        return primitive(js, tok, mx);
    } else if (tok->type == JSMN_STRING) {
        return value(js, tok, mx);
    } else if (tok->type == JSMN_OBJECT) {
        return object(js, tok, mx);
    } else if (tok->type == JSMN_ARRAY) {
        return array(js, tok, mx);
    }
    return 0;
}

static jsmntok_t * parse(const char *js, size_t jslen) {
    int r;
    jsmn_parser p;
    jsmntok_t *tok = NULL;
    size_t tokcount = 2;
    
    jsmn_init(&p);
    tok = mxMalloc(sizeof(*tok) * tokcount);
    if (tok == NULL) {
        mexErrMsgTxt("mxMalloc()");
    }
    
    for (;;) {
        r = jsmn_parse(&p, js, jslen, tok, tokcount);
        if (r < 0) {
            if (r == JSMN_ERROR_NOMEM) {
                tokcount = tokcount * 2;
                tok = mxRealloc(tok, sizeof(*tok) * tokcount);
                if (tok == NULL) {
                    mexErrMsgTxt("mxRealloc()");
                }
            }
            else if ((r == JSMN_ERROR_INVAL) || (r == JSMN_ERROR_PART)) {
                mexErrMsgTxt("Invalid or incomplete JSON.");
            }
            else {
                mexErrMsgTxt("Unknown JSON parsing error.");
            }
        }
        else {
            break;
        }
    }
    
    return tok;
}

static char * get_data(const mxArray * mx, size_t * jslen) {
    /* should attempt to minimise copy */
    int i, filename, sts;
    mxArray *ma = NULL;
    char *js = NULL;

    js = mxArrayToString(mx);
    if (js == NULL) {
        mexErrMsgTxt("mxArrayToString()");
    }
    *jslen = strlen(js);
    if (*jslen == 0)
        mexErrMsgTxt("Empty JSON.");
    
    /* detect whether input string is a filename */
    for (i = 0, filename = 1; i < *jslen; i++) {
        if ((js[i] == '{') || (js[i] == '[')) {
            filename = 0;
            break;
        }
    }
    if (filename == 1) {
        mxFree(js);
        sts = mexCallMATLAB(1, &ma, 1, (mxArray **)&mx, "fileread");
        if (sts != 0) {
            mexErrMsgTxt("Cannot read JSON file.");
        }
        js = mxArrayToString(ma);
        if (js == NULL) {
            mexErrMsgTxt("mxArrayToString()");
        }
        mxDestroyArray(ma);
        *jslen = strlen(js);
    }
    return js;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    char *js = NULL;
    const char *field = NULL;
    size_t i, jslen = 0, nfields;
    jsmntok_t *tok = NULL;
    mxArray *mx = NULL;
    char *repsty = NULL;

    /* Validate input arguments */
    if (nrhs == 0) {
        mexErrMsgTxt("Not enough input arguments.");
    }
    else if (nrhs > 2) {
        mexErrMsgTxt("Too many input arguments.");
    }
    if (!mxIsChar(prhs[0])) {
        mexErrMsgTxt("Input must be a string.");
    }
    ReplacementStyle = JSON_REPLACEMENT_STYLE_UNDERSCORE;
    if (nrhs > 1) {
        if (!mxIsStruct(prhs[1])){
            mexErrMsgTxt("Input must be a struct.");
        }
        nfields = mxGetNumberOfFields(prhs[1]);
        for (i = 0; i < nfields; i++) {
            field = mxGetFieldNameByNumber(prhs[1], i);
            if (!strcasecmp(field,"replacementStyle")) {
                mx = mxGetFieldByNumber(prhs[1],0,i);
                if (mx != NULL) {
                    repsty = mxArrayToString(mx);
                    if (!strcasecmp(repsty,"nop")) {
                        ReplacementStyle = JSON_REPLACEMENT_STYLE_NOP;
                    }
                    else if (!strcasecmp(repsty,"underscore")) {
                        ReplacementStyle = JSON_REPLACEMENT_STYLE_UNDERSCORE;
                    }
                    else if (!strcasecmp(repsty,"hex")) {
                        ReplacementStyle = JSON_REPLACEMENT_STYLE_HEX;
                    }
                    else if (!strcasecmp(repsty,"delete")) {
                        ReplacementStyle = JSON_REPLACEMENT_STYLE_DELETE;
                    }
                    else {
                        mexErrMsgTxt("Unknown replacementStyle.");
                    }
                    mxFree(repsty);
                }
            }
            else {
                mexErrMsgTxt("Unknown optional parameter.");
            }
        }
    }

    /* Get JSON data as char array */
    js = get_data(prhs[0], &jslen);

    /* Parse JSON data */
    tok = parse(js, jslen);

    /* Create output structure */
    create_struct(js, tok, &plhs[0]);

    mxFree(js);
    mxFree(tok);
}
