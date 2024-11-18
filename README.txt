This is the SSPsyGene knowledge base demo / testing code.

Not hard to install:
- No dependencies
- Copy to any cgi-bin directory of a webserver and run it. 
- In Apache: directory should have these lines:
    Options Indexes FollowSymLinks ExecCGI MultiViews
    AddHandler cgi-script .cgi .py
    AllowOverride None
    Require all granted

Then open show.py in your internet browser.


Data wrangling is under data/ 
See README and log.txt files there
