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
      <div class="col-xs-8" role="main">
	
      	<table class="table table-striped table-condensed table-hover">
	  <thead><tr>
	      <th><code>{{filename}}</code></th>
	      <th><button type="button" class="btn btn-block btn-xs btn-primary">
		  Platform</button></th>
	      <th><button type="button" class="btn btn-block btn-xs btn-primary">
		  Test</button></th>
	      <th><button type="button" class="btn btn-block btn-xs btn-primary">
		  Result</button></th>
	  </tr></thead>

	  <tbody>
	    % for res in testres:
	    %  testcode=res['resultcode']
	    <tr>
	      <td></td>
	      <td>{{res['platform']}}</td>
	      <td>{{res['test']}}</td>
	      <td><a href="{{res['url']}}">
		% if testcode == 0:
		  <button class="btn btn-block btn-xs btn-success">
		% elif testcode == 1:
		  <button class="btn btn-block btn-xs btn-info">
		% elif testcode == 2:
		  <button class="btn btn-block btn-xs btn-warning">
		% elif testcode == 3:
		  <button class="btn btn-block btn-xs btn-danger">
		% end
		{{res['result']}}</button>
	      </a></td>
	    </tr>
	    % end
	  </tbody>
	</table>
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
