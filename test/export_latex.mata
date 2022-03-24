cap mata mata drop export_latex()
mata
void function export_latex(
    real matrix M,
    string scalar out,
    string vector header,
    string vector footer,
    string vector collabels,
    string vector rowlabels,
    string matrix rowformats,
    real   vector lindex,
    string vector lines)
{
    string scalar row
    string colvector Mfmt

    real scalar i, j, fh
    real scalar labrow, labcol, fmtrow, fmtcol, head, foot

    labrow = any(rowlabels  :!= "")
    labcol = any(collabels  :!= "")
    fmtrow = any(rowformats :!= "")
    fmtcol = 0
    head   = any(header :!= "")
    foot   = any(footer :!= "")

    if (labrow) {
        if (rows(M) != length(rowlabels)) {
            errprintf("Matrix had %g rows and labels had %g\n", cols(M), length(rowlabels))
            error(198)
        }
    }

    if (fmtrow) {
        if (rows(M) == length(rowformats)) {
            fmtrow = 1
            fmtcol = 0
        }
        else if (cols(M) == length(rowformats)) {
            fmtrow = 0
            fmtcol = 1
            printf("Formats interpreted as column formats\n")
        }
        else if ( (rows(M) == rows(rowformats)) & (cols(M) == cols(rowformats)) ) {
            fmtrow = 1
            fmtcol = 1
            printf("Formats interpreted as individual cell formats\n")
        }
        else {
            errprintf("Matrix had %g rows and formats had %g\n", cols(M), length(rowformats))
            error(198)
        }
    }

    if (labcol) {
        if (cols(M) != length(collabels)) {
            errprintf("Matrix had %g columns and labels had %g\n", cols(M), length(collabels))
            error(198)
        }

        row = labrow? "&": ""
        for (j = 1; j < cols(M); j++) {
            row = row + collabels[j] + "&"
        }
        row = row + collabels[j] + "\\"
    }

    Mfmt = labcol? row: J(0, 1, "")

    for(i = 1; i <= rows(M); i++) {
        row = labrow? "                " + rowlabels[i] + "&": ""
        for(j = 1; j < cols(M); j++) {
            if ( missing(M[i, j]) ) {
                row = row + "&"
            }
            else {
                if ( (fmtrow == 1) & (fmtcol == 0) ) {
                    row = row + sprintf(rowformats[i], M[i, j]) + "&"
                }
                else if ( (fmtrow == 0) & (fmtcol == 1) ) {
                    row = row + sprintf(rowformats[j], M[i, j]) + "&"
                }
                else if ( (fmtrow == 1) & (fmtcol == 1) ) {
                    row = row + sprintf(rowformats[i, j], M[i, j]) + "&"
                }
                else {
                    row = row + sprintf("%12.0g", M[i, j]) + "&"
                }
            }
        }

        if ( missing(M[i, j]) ) {
            row = row + "\\"
        }
        else {
            if ( (fmtrow == 1) & (fmtcol == 0) ) {
                row = row + sprintf(rowformats[i], M[i, j]) + "\\"
            }
            else if ( (fmtrow == 0) & (fmtcol == 1) ) {
                row = row + sprintf(rowformats[j], M[i, j]) + "\\"
            }
            else if ( (fmtrow == 1) & (fmtcol == 1) ) {
                row = row + sprintf(rowformats[i, j], M[i, j]) + "\\"
            }
            else {
                row = row + sprintf("%12.0g", M[i, j]) + "\\"
            }
        }
        Mfmt = Mfmt \ row
    }

    for(i = 1; i <= length(lindex); i++) {
        Mfmt[lindex[i]] = Mfmt[lindex[i]] + lines[i]
    }

    if (head) {
        Mfmt = header \ Mfmt
    }

    if (foot) {
        Mfmt = Mfmt \ footer
    }

    fh = fopen(out, "rw")
    for(i = 1; i <= rows(Mfmt); i++) {
        fwrite(fh, Mfmt[i])
        fwrite(fh, sprintf("\n"))
    }
    fclose(fh)
}
end
