<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">

<head>
	<meta http-equiv="content-type" content="text/html; charset=iso-8859-1" />
	<title>Upload a presentation file</title>
</head>

<body>

<?php  //script 11.4 in Ullman PHP text:  upload_file.php
// this page uploads an arbitrary file to a directory

// address error handling
ini_set ('display_errors', 1); // turn on error display
error_reporting (E_ALL & ~E_NOTICE);

if (isset ($_POST['submit'])) { // handle form

	// try to move the uploaded file

	if (move_uploaded_file ($_FILES['thefile']['tmp_name'], "Uploads/{$_FILES['thefile']['name']}")) {

		print '<p>Your file has been uploaded.</p>';

	} else { // error handling 

		print '<p>Your file could not be uploaded because:<b>';

		// Print error message
		switch ($_FILES['thefile']['error']) {
			case 1:
				print 'The file exceeds the upload_max_filesize setting in php.ini.';
				break;
			case 2:
				print 'The file exceeds the MAX_FILE_SIZE setting in the html form.';
				break;
			case 3:
				print 'The file was only partially uploaded.';
				break;
			case 4:
				print 'No file was uploaded.';
				break;
		}
		print '</b>.</p>';
	}
} // End of submit if


# end of php
?>

<form action="submit_hw_file.php" enctype="multipart/form-data" method="post"> 
<p>You may upload your presentation file using this form.  
<br/>
The file must be in pdf format. 
It should have a filename that is of the form: uinetid_ps#_q#.pdf, where uinetid
is your official uinetid in lower case letters, and  the #'s should be replaced by
the number of the problem set, and the question your were assigned respectively.
For PS3, the question # should be odd, that is those assign 7-8 should have "q7".
You are free to use any text formatting system you want to prepare the pdf
file, but eventually you will probably want to learn latex and for latex it
is quite convenient to use beamer.  You can google to find useful tutorials.
<br/>
Files that are larger than 3Mbytes will be rejected.  
<br/><br />

<input type="hidden" name="MAX_FILE_SIZE" value="3000000" />
<input type="file" name="thefile" /><br /><br />
<input type="submit" name="submit" value="Upload This File" />
</p>
</form>
</body>
</html>
