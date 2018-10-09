#!/usr/bin/env bash

for Name in EXSdetect_Download EXSdetect_Input EXSdetect_Output EXSdetect_Tricks SweepLine_Introduction SweepLine_Basic SweepLine_More FOF_Introduction
do
	sed -n '1,/<!--ENDHEAD-->/p' EXSdetect_Introduction.html |\
		sed "/${Name}/s/<\\!--current_page_item-->/current_page_item/" >${Name}.new.html
	sed -n '/<!--BEGINCONTENT-->/,/<!--ENDCONTENT-->/p' ${Name}.html >> ${Name}.new.html
	sed -n '/<!--BEGINTAIL-->/,$p' EXSdetect_Introduction.html >> ${Name}.new.html

	\mv ${Name}.new.html ${Name}.html
done
