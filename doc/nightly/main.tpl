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
	<option value="a">File type</option>
	% for t in types:
	    <option value="{{t}}">{{t}}</option>
	% end
    </select></th>
    <th><select id="version">
	<option value="a">Version</option>
	% for v in versions:
	    <option value="{{v}}">{{v}}</option>
	% end
    </select></th>
    <th><select id="branch">
	<option value="all">Branch</option>
	% for b in branches:
	    <option value="{{b}}">{{b}}</option>
	% end
    </select></th>
    <th><button id="commit">Commit</button></th>
    <th><button id="date">Uploaded</button></th>
    <th><button id="size">File size</button></th>
  </tr></thead>
  <tbody>
    % for file in files:
    %   dir=urlmap[file['type']]
    %   hash=file['hash']
        <tr>
	  <td><input type="button" value="&#x25BC;" 
		     onclick="location.href='/get/{{dir}}/{{hash}}';"/></td>
	  <td>{{file['type']}}</td>
	  <td>{{file['version']}}</td>
	  <td>{{file['branch']}}</td>
	  <td><code>{{file['hash']}}</code></td>
	  <td>{{file['date']}}</td>
	  <td>{{file['size']}}</td>
       </tr>
    % end
  </tbody>
</table>
</form>
