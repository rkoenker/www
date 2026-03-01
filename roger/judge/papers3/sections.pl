# LaTeX2HTML 96.1 (Feb 5, 1996)
# Associate sections original text with physical files.

$key = q/0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0/;
$section_info{$key} = '0%:%papers3.html%:%No Title' unless ($section_info{$key}); 
$done{"papers3.html"} = 1;
$key = q/0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0/;
$section_info{$key} = '3%:%node1.html%:%  About this document ... ' unless ($section_info{$key}); 
$done{"node1.html"} = 1;

1;

