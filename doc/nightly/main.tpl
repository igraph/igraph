<!doctype html> <!-- -*- html-mode -*- -->
<html>
<head>
  <link rel="stylesheet" href="style.css" type="text/css" />
  <title>igraph nightly builds</title>
</head>

<body>

<h2>igraph nightly builds</h2>

<h3>Latest files</h3>

<p>Select the file type and also the branch you want. If you want 
	a snapshot of the next major release, then choose
	the development branch.
</p>

<form>
<table><thead><tr>
    <th><!-- download buttons --></th>
    <th><select id="type">
	<option value="all">File type (all)</option>
	% for t in types:
	    <option value="{{t}}">Type {{urlmap[t]}}</option>
	% end
    </select></th>
    <th><select id="version">
	<option value="all">Version (all)</option>
	% for v in versions:
	    <option value="{{v}}">Version {{v}}</option>
	% end
    </select></th>
    <th><select id="branch">
	<option value="all">Branch (all)</option>
	% for b in branches:
	    <option value="{{b}}">Branch {{b}}</option>
	% end
    </select></th>
    <th><button id="commit">Commit</button></th>
    <th><button id="date">Uploaded</button></th>
    <th><button id="size">File size</button></th>
  </tr></thead>
  <tbody>
    % for file in files:
    %   filename=file['filename']
        <tr>
	  <td><input type="button" value="&#x25BC;" 
		     onclick="location.href='/get/{{filename}}';"/></td>
	  <td>{{urlmap[file['type']]}}</td>
	  <td>{{file['version']}}</td>
	  <td>{{file['branch']}}</td>
	  <td><code>{{file['hash']}}</code></td>
	  <td>{{file['date']}}</td>
	  <td>{{file['size']}}</td>
       </tr>
    % end
  </tbody>
</table>
<noscript><input type="submit" value="Submit"></noscript>
</form>

<script>

document.getElementById('type').value="{{dtype}}";
document.getElementById('version').value="{{version}}";
document.getElementById('branch').value="{{branch}}";

function redirect() {
  var type=document.getElementById('type').value || "all";
  var version=document.getElementById('version').value || "all";
  var branch=document.getElementById('branch').value || "all";
  var goto = "/list/" + type + "/" + version + "/" + branch
  window.location = goto;
};
document.getElementById('type').onchange = redirect;
document.getElementById('version').onchange = redirect;
document.getElementById('branch').onchange = redirect;

</script>
