#!/bin/bash
#SBATCH --job-name=getfasta
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH -p debug
##SBATCH --qos=fast
#SBATCH --cpus-per-task=4
#SBATCH --array=1-94


#################################################################################

OUTPUT='/home/newatil/fasta_files_merged'
mkdir -p ${OUTPUT}
cd ${OUTPUT}

##############################################################################

#genes=('FPM_i_1_A' 'FPM_i_2_A' 'FPM_m_1_D' 'FPM_m_3_A' 'FPM_m_3_D' 'FPM_m_4_B' 'FPM_m_4_D' 'FPM_m_5_B' 'FPM_m_5_D' 'FPM_m_6_D' 'FPM_m_7_D' 'FPM_m_8_D' 'HPM_1Ax' 'HPM_1Ay' 'HPM_1Bx' 'HPM_1By' 'HPM_1Dx' 'HPM_1Dy')

#var=('bc2001' 'bc2002' 'bc2003' 'bc2004' 'bc2005' 'bc2006' 'bc2007' 'bc2008' 'bc2009' 'bc2010' 'bc2011' 'bc2012' 'bc2013' 'bc2014' 'bc2015' 'bc2016' 'bc2017' 'bc2018' 'bc2019' 'bc2020' 'bc2021' 'bc2022' 'bc2023' 'bc2024' 'bc2025' 'bc2026' 'bc2027' 'bc2028' 'bc2029' 'bc2030' 'bc2032' 'bc2033' 'bc2034' 'bc2035' 'bc2036' 'bc2037' 'bc2038' 'bc2039' 'bc2040' 'bc2041' 'bc2042' 'bc2043' 'bc2044' 'bc2045' 'bc2046' 'bc2047' 'bc2048' 'bc2049' 'bc2050' 'bc2051' 'bc2052' 'bc2053' 'bc2054' 'bc2055' 'bc2056' 'bc2057' 'bc2058' 'bc2059' 'bc2060' 'bc2061' 'bc2062' 'bc2063' 'bc2064' 'bc2065' 'bc2066' 'bc2067' 'bc2068' 'bc2069' 'bc2070' 'bc2071' 'bc2072' 'bc2073' 'bc2074' 'bc2075' 'bc2076' 'bc2077' 'bc2078' 'bc2079' 'bc2080' 'bc2082' 'bc2083' 'bc2084' 'bc2085' 'bc2086' 'bc2087' 'bc2088' 'bc2089' 'bc2090' 'bc2091' 'bc2092' 'bc2093' 'bc2094' 'bc2095' 'bc2096')


echo '**************************** PROGRAM START RUNING ******************************************'

for v in 'bc2001' 'bc2002' 'bc2003' 'bc2004' 'bc2005' 'bc2006' 'bc2007' 'bc2008' 'bc2009' 'bc2010' 'bc2011' 'bc2012' 'bc2013' 'bc2014' 'bc2015' 'bc2016' 'bc2017' 'bc2018' 'bc2019' 'bc2020' 'bc2021' 'bc2022' 'bc2023' 'bc2024' 'bc2025' 'bc2026' 'bc2027' 'bc2028' 'bc2029' 'bc2030' 'bc2032' 'bc2033' 'bc2034' 'bc2035' 'bc2036' 'bc2037' 'bc2038' 'bc2039' 'bc2040' 'bc2041' 'bc2042' 'bc2043' 'bc2044' 'bc2045' 'bc2046' 'bc2047' 'bc2048' 'bc2049' 'bc2050' 'bc2051' 'bc2052' 'bc2053' 'bc2054' 'bc2055' 'bc2056' 'bc2057' 'bc2058' 'bc2059' 'bc2060' 'bc2061' 'bc2062' 'bc2063' 'bc2064' 'bc2065' 'bc2066' 'bc2067' 'bc2068' 'bc2069' 'bc2070' 'bc2071' 'bc2072' 'bc2073' 'bc2074' 'bc2075' 'bc2076' 'bc2077' 'bc2078' 'bc2079' 'bc2080' 'bc2082' 'bc2083' 'bc2084' 'bc2085' 'bc2086' 'bc2087' 'bc2088' 'bc2089' 'bc2090' 'bc2091' 'bc2092' 'bc2093' 'bc2094' 'bc2095' 'bc2096'; do for g in 'FPM_i_1_A' 'FPM_i_2_A' 'FPM_m_1_D' 'FPM_m_3_A' 'FPM_m_3_D' 'FPM_m_4_B' 'FPM_m_4_D' 'FPM_m_5_B' 'FPM_m_5_D' 'FPM_m_6_D' 'FPM_m_7_D' 'FPM_m_8_D' 'HPM_1Ax' 'HPM_1Ay' 'HPM_1Bx' 'HPM_1By' 'HPM_1Dx' 'HPM_1Dy';do cat demux.${v}_${g}_test_stplus.fasta  demux.${v}_${g}_test_stminus.rc.fasta > demux.${v}_${g}_test_merged.fasta;done;done

echo '************************** DONE ***********************************************'
