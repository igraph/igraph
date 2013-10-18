<!doctype html> <!-- -*- html-mode -*- -->
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">

  <link href="http://www.igraph.org/bootstrap/css/bootstrap.css"
	rel="stylesheet">
  <link href="http://www.igraph.org/css/other.css"
	rel="stylesheet">
  <link href='http://fonts.googleapis.com/css?family=Lato:300,400,700'
        rel='stylesheet' type='text/css'>
  <link href="//netdna.bootstrapcdn.com/font-awesome/3.2.1/css/font-awesome.css"
	rel="stylesheet">

  <link rel="stylesheet" href="style.css" type="text/css" />

  <title>igraph nightly builds</title>
</head>

<body>

  <body>

    <div>
      <a href="https://github.com/igraph/igraph">
        <img style="position: absolute; top:50px; right: 0; border: 0; z-index:1000" 
             src="https://s3.amazonaws.com/github/ribbons/forkme_right_red_aa0000.png"
             alt="Fork me on GitHub">
      </a>
    </div>

    <div class="navbar navbar-default navbar-fixed-top">
      <div class="container">
        <div class="navbar-header">
          <button type="button" class="navbar-toggle"
                  data-toggle="collapse" data-target=".navbar-collapse">
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a class="navbar-brand" href="#">
            <img src="http://www.igraph.org/img/igraph2.svg" width=25>
            igraph
          </a>
        </div>
        <div class="navbar-collapse collapse">
          <ul class="nav navbar-nav navbar-left">
            <li class="active"><a href="#">Home</a></li>
            <li><a href="">News</a></li>
            <li><a href="">Get started</a></li>
            <li><a href="">Download</a></li>
            <li><a href="">Get help</a></li>
            <li><a href="">Contribute</a></li>
          </ul>
          <form class="navbar-form navbar-right">
            <div class="form-group">
              <input type="text" placeholder="Search documentation"
                     class="form-control" size="25">
            </div>
            <button type="submit" class="btn btn-success">
              <span class="glyphicon glyphicon-search"></span></button>
          </form>
        </div><!--/.navbar-collapse -->
      </div>
    </div>

    <div class="jumbotron">
      <div class="container">
        <h1 class="hidden-xs mainheader">
          <img src="http://www.igraph.org/img/igraph2.svg"> igraph
          nightly builds
        </h1>
        <h1 class="visible-xs">
            igraph<br/>
	    nightly builds
        </h1>
        <p class="lead">
	  Nightly builds from the branches under active development
        </p>
      </div>
    </div>

    <div class="container">
      <div class="col-xs-12 bs-docs-section" role="main">

	<form>
	  <table class="table table-striped table-condensed"><thead><tr>
		<th><!-- download buttons --></th>
		<th><select id="type" class="form-control btn-primary">
		    <option value="all">File type (all)</option>
		    % for t in types:
		    <option value="{{t}}">Type {{urlmap[t]}}</option>
		    % end
		</select></th>
		<th><select id="version" class="form-control btn-primary">
		    <option value="all">Version (all)</option>
		    % for v in versions:
		    <option value="{{v}}">Version {{v}}</option>
		    % end
		</select></th>
		<th><select id="branch" class="form-control btn-primary">
		    <option value="all">Branch (all)</option>
		    % for b in branches:
		    <option value="{{b}}">Branch {{b}}</option>
		    % end
		</select></th>
		<th><button id="commit" class="btn btn-primary">
		    Commit</button></th>
		<th><button id="date" class="btn btn-primary">
		    Uploaded</button></th>
		<th><button id="size" class="btn btn-primary">
		    File size</button></th>
	    </tr></thead>
	    <tbody>
	      % for file in files:
	      %   filename=file['filename']
              <tr>
		<td><a href="/get/{{filename}}">
		    <i class="icon-download"></i></a></td>
		<td>{{urlmap[file['type']]}}</td>
		<td>{{file['version']}}</td>
		<td><a href="https://github.com/igraph/igraph/tree/{{file['branch']}}">{{file['branch']}}</a></td>
		<td><a href="https://github.com/igraph/igraph/commits/{{file['hash']}}"><code>{{file['hash']}}</code></a></td>
		<td>{{file['date']}}</td>
		<td>{{human_size(file['size'])}}</td>
	      </tr>
	      % end
	    </tbody>
	  </table>
	  <noscript><input type="submit" value="Submit"></noscript>
	</form>
      </div>
    </div>

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

    <div class="container footer">
      <div class="row">
        <div class="col-xs-6">
        </div>
        <div class="col-xs-6 text-right">
          &copy; 2013 The igraph core team
        </div>
      </div>
    </div>

    <!-- Bootstrap core JavaScript
    ================================================== -->
    <!-- Placed at the end of the document so the pages load faster -->
    <script src="js/jquery-2.0.3.min.js"></script>
    <script src="bootstrap/js/bootstrap.min.js"></script>
    <script src="js/affix.js"></script>
    <script type="text/javascript">
      $(document).ready(function(){
        $( ".thumbnail" ).mouseenter(function() {
          $(this).find('.caption').removeClass("flipOutY").addClass("flipInY").show();
        })
        .mouseleave(function() {
          $(this).find('.caption').removeClass("flipInY").addClass("flipOutY");
});
})
    </script>
  </body>
</html>
