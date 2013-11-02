% from bottle import template
% include('header.tpl')
% include('github.tpl')
% include('navbar.tpl')
% include('ijumbo.tpl')

    <div class="container">
      <div class="col-xs-12" role="main">

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
		      <li><a href="/{{url}}/all/{{version}}/{{branch}}">
			  All types
			  % if dtype=="all":
			    <i class="icon-ok"></i>
			  % end
		      </a></li>
		      <li class="divider"></li>
		      % for t in types:
		        <li><a href="/{{url}}/{{t}}/{{version}}/{{branch}}">
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
		      <li><a href="/{{url}}/{{dtype}}/all/{{branch}}">
			  All versions
			  % if version=="all":
			    <i class="icon-ok"></i>
			  % end
		      </a></li>
		      <li class="divider"></li>
		      % for v in versions:
		        <li><a href="/{{url}}/{{dtype}}/{{v}}/{{branch}}">
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
		      <li><a href="/{{url}}/{{dtype}}/{{version}}/all">
			  All branches
			  % if branch=="all":
			    <i class="icon-ok"></i>
			  % end
		      </a></li>
		      <li class="divider"></li>
		      % for b in branches:
		        <li><a href="/{{url}}/{{dtype}}/{{version}}/{{b}}">
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
      </div>
    </div>

% include('footer.tpl')
