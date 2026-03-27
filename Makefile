notebook:
	claude --system-prompt-file agents/author.prompt --tools "" -p "Generate the notebook as specified." > sandbox/monomer_filter_cluster.py
	cat sandbox/monomer_filter_cluster.py | claude --system-prompt-file agents/reviewer.prompt --tools "" -p "Review and correct the above marimo notebook." > /tmp/reviewed_notebook.py
	mv /tmp/reviewed_notebook.py sandbox/monomer_filter_cluster.py