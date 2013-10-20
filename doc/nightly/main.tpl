<!doctype html> <!-- -*- html-mode -*- -->
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">

  <link href="http://www.igraph.org/bootstrap/css/bootstrap.css"
	rel="stylesheet">
  <link href='http://fonts.googleapis.com/css?family=Lato:300,400,700'
        rel='stylesheet' type='text/css'>
  <link href="http://netdna.bootstrapcdn.com/font-awesome/3.2.1/css/font-awesome.css"
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
          <a class="navbar-brand" href="http://www.igraph.org">
            <img src="http://www.igraph.org/img/igraph2.svg" width=25>
            igraph
          </a>
        </div>
        <div class="navbar-collapse collapse">
          <ul class="nav navbar-nav navbar-left">
            <li><a href="http://www.igraph.org">Home</a></li>
            <li><a href="http://www.igraph.org/news">News</a></li>
            <li><a href="http://www.igraph.org/all.html">Get started</a></li>
            <li><a href="http://www.igraph.org/all.html#downloads">
		Download</a></li>
            <li><a href="http://www.igraph.org/all.html#help">Get help</a></li>
            <li><a href="http://www.igraph.org/all.html#contribute">
		Contribute</a></li>
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
	  Nightly builds from branches under active development.
        </p>
      </div>
    </div>

    <div class="container">
      <div class="col-xs-12" role="main">

	<form>
	  <table class="table table-striped table-condensed table-hover">
	    <thead><tr>
		<th><!-- download buttons --></th>

		<th><div class="btn-group btn-block">
		    <button type="button"
			    class="btn btn-primary btn-block dropdown-toggle"
			    data-toggle="dropdown">File type
		      % if dtype != "all":
		        <i class="icon-filter"></i>
		      % end
		    </button>
		    <ul class="dropdown-menu" role="menu">
		      <li><a href="/list/all/{{version}}/{{branch}}">
			  All types
			  % if dtype=="all":
			    <i class="icon-ok"></i>
			  % end
		      </a></li>
		      <li class="divider"></li>
		      % for t in types:
		        <li><a href="/list/{{t}}/{{version}}/{{branch}}">
			    {{urlmap[t]}}
			% if dtype==t:
			    <i class="icon-ok"></i>
			% end
			</a></li>
		      % end
		    </ul>
		</div></th>

		<th><div class="btn-group btn-block">
		    <button type="button"
			    class="btn btn-primary btn-block dropdown-toggle"
			    data-toggle="dropdown">Version
		      % if version != "all":
		        <i class="icon-filter"></i>
		      % end
		    </button>
		    <ul class="dropdown-menu" role="menu">
		      <li><a href="/list/{{dtype}}/all/{{branch}}">
			  All versions
			  % if version=="all":
			    <i class="icon-ok"></i>
			  % end
		      </a></li>
		      <li class="divider"></li>
		      % for v in versions:
		        <li><a href="/list/{{dtype}}/{{v}}/{{branch}}">
			    {{v}}
			% if version==v:
			    <i class="icon-ok"></i>
			% end
			</a></li>
		      % end
		  </ul>
		</div></th>

		<th><div class="btn-group btn-block">
		    <button type="button"
			    class="btn btn-primary btn-block dropdown-toggle"
			    data-toggle="dropdown">Branch
		      % if branch != "all":
		        <i class="icon-filter"></i>
		      % end
		    </button>
		    <ul class="dropdown-menu" role="menu">
		      <li><a href="/list/{{dtype}}/{{version}}/all">
			  All branches
			  % if branch=="all":
			    <i class="icon-ok"></i>
			  % end
		      </a></li>
		      <li class="divider"></li>
		      % for b in branches:
		        <li><a href="/list/{{dtype}}/{{version}}/{{b}}">
			    {{b}}
			% if branch==b:
			    <i class="icon-ok"></i>
			% end
			</a></li>
		      % end
		    </ul>
		</div></th>

		<th><button type="button" class="btn btn-block btn-primary">
		    Commit</button></th>
		<th><button type="button" class="btn btn-block btn-primary">
		    Uploaded</button></th>
		<th><button type="button" class="btn btn-block btn-primary">
		    File size</button></th>
		<th><button type="button" class="btn btn-block btn-primary">
		    Tests</button></th>

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
		<td><a href="/tests/{{filename}}">
		  % if filename in tests:
		  %   testcode=long(tests[filename])
		  %   if testcode == 0:
		        <button class="btn btn-block btn-xs btn-success">
			  OK</button>
		  %   elif testcode == 1:
			<button class="btn btn-block btn-xs btn-info">
			  NOTE</button>
		  %   elif testcode == 2:
			<button class="btn btn-block btn-xs btn-warning">
			  WARNING</button>
		  %   elif testcode == 3:
			<button class="btn btn-block btn-xs btn-danger">
			  ERROR</button>
                  %   end
		  % else:
			<!-- Untested -->
		  % end
		</a></td>
	      </tr>
	      % end
	    </tbody>
	  </table>
	  <noscript><input type="submit" value="Submit"></noscript>
	</form>
      </div>
    </div>

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
    <script src="http://www.igraph.org/js/jquery-2.0.3.min.js"></script>
    <script src="http://www.igraph.org/bootstrap/js/bootstrap.min.js"></script>
    <script>
      $('.dropdown-toggle').dropdown()
    </script>
  </body>
</html>
