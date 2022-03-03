cap program drop output
program output
	syntax, treatment(varname)                     /// control variable
            MATAsave(str)                       /// input saved mata object
     [                                          ///
		outpath(str)							/// output filepath/name
     ]

	glevelsof `treatment', loc(treatlvls)
	local treatarms = `:list sizeof treatlvls' - 1

	mata panelA_matrix = J(`= `treatarms' * 2', `=(`treatarms' - 1)*3 + 2', .)

	mata panelB_matrix = J(`= `treatarms' * 3', `=`treatarms' + 1', .)

	forval i=1/`treatarms' {
		mata panelA_matrix[`=(`i' - 1)*2 + 1', .] = `matasave'.decomposition.est[`i', .]
		mata panelA_matrix[`=(`i' - 1)*2 + 2', .] = `matasave'.decomposition.se[`i', .]
		mata panelB_matrix[`=(`i' - 1)*3 + 1', .] = `matasave'.estimates.est[`i', .]
		mata panelB_matrix[`=(`i' - 1)*3 + 2', .] = `matasave'.estimates.se_po[`i', .]
		mata panelB_matrix[`=(`i' - 1)*3 + 3', .] = `matasave'.estimates.se_or[`i', .]
	}

	mata ncolsA = cols(panelA_matrix)
	mata ncolsB = cols(panelB_matrix)
	mata temp = (ncolsA - ncolsB) / 2
	mata blank = J(rows(panelB_matrix), temp, .)
	mata panelB_matrix = blank, panelB_matrix, blank

	* 2. make rowlabels mata object and other helper objects
	mata table1_matrix = panelA_matrix \ panelB_matrix
	
	mata panelA_rlabels = /// 
		"Small Class Size" \ ///
		"" \ ///
		"Teaching Aide" \ ///
		""
		
	mata panelB_rlabels = /// 
		"Small Class Size" \ ///
		"" \ ///
		"" \ ///
		"Teaching Aide" \ ///
		"" \ ///
		"" 
		
	mata table1_rlabels = panelA_rlabels \ panelB_rlabels
	
	mata panelA_rfmt = ///
		"%5.3f" \ ///
		"(%5.3f)" \ ///
		"%5.3f" \ ///
		"(%5.3f)"
	
	mata panelB_rfmt = ///
		"%5.3f" \ ///
		"(%5.3f)" \ ///
		"[%5.3f]" \ ///
		"%5.3f" \ ///
		"(%5.3f)" \ ///
		"[%5.3f]"
		
	mata table1_rfmt = panelA_rfmt \ panelB_rfmt
	
	mata panelA_header = ///
		"\begin{table}[h]" \ ///
		"\caption{\label{tab:contamination}Project STAR Contamination Bias and Treatment Effect Estimates}" \ ///
		"\vspace{-0.25cm}" \ ///
		"\begin{center}" \ ///
		"\addtolength{\leftskip} {-2cm}" \ ///
		"\addtolength{\rightskip}{-2cm}" \ ///
		"\begin{tabular}{l c>{\hspace*{-1mm}}c>{\hspace*{-3.5mm}}c>{\hspace*{-2mm}}c>{\hspace*{1mm}}c}" \ ///
		"\toprule " \ ///
		"&\multicolumn{5}{c}{A. Contamination Bias Estimates}\\" \ ///
		"\cline{2-6}&Regression&Own&\multirow{2}{*}{Bias}&\multicolumn{2}{c}{Worst-Case Bias}\\" \ ///
		"\cline{5-6}&Coefficient&Effect&&Negative&Positive\\" \ ///
		"&(1)&(2)&(3)&(4)&(5)\\" \ ///
		"\midrule"
		
	mata panelB_header = ///
		"\midrule" \ ///
		"&\multicolumn{5}{c}{B. Treatment Effect Estimates}\\" \ ///
		"\cline{2-6}&&Unweighted&\multicolumn{2}{c}{Efficiently-Weighted}&\\" \ ///
		"\cline{4-5}&&(ATE)&One-at-a-time&Common&\\" \ ///
		"&&(1)&(2)&(3)&\\" \ ///
		"\midrule"
		
	mata panelA_footer = ///
		" " \ ///
		" "
	
	mata panelB_footer = ///
		"\bottomrule" \ ///
		"\end{tabular}" \ ///
		"\par\end{center}" \ ///
		"{\footnotesize{}}" \ ///
		"{\footnotesize\par}" \ ///
		"\end{table}" \ ///
		"" \ ///
		"" \ ///
		"" \ ///
		""
	
	mata panelB_insert = 4 \ 4 \ 4 \ 4 \ 4 \ 4
	
	* 3. use export_latex() mata function to export latex table!
	if (`"`outpath'"' == "") local outpath table1
	local out "`outpath'.tex"
	
	mata export_latex(table1_matrix, `"`out'"', panelA_header, panelB_footer, "", table1_rlabels, table1_rfmt, panelB_insert, panelB_header)
	
end
