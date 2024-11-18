#!/usr/bin/env python3

# this can be run from the command line with CGI arguments:
# /usr/bin/python3 show.py 'filter_tom=on&filter_syndromic=on'

# format python errors in html, as we're a CGI script
import cgi, timeit, sys, os, sqlite3
import cgitb; cgitb.enable()

## these are default python modules, no errors expected here
#import urllib, urllib2, zlib, collections, StringIO, gzip, string, sys, os, random, \
    #subprocess, re, types, socket, cPickle, copy, glob, tempfile, Cookie

from collections import defaultdict, namedtuple
from os.path import *
from distutils import spawn

DEBUG = False
cgiArgs = None

#tableEngine = "tsv"
#tableInput = "genes.tsv"
#tableEngine = "sqlite"

# short table descriptions
tableShort = {
    "mousePerturb4Tf" : "MouseGene expr",
    "comp" : "Mouse Cell Type Comp",
    "zebraAsdSizes" : "Zebra Brain Sizes",
    "perturbFishAstro" : "Perturb-Fish Gene Expr Astrocytes"
}

tableDescs = {
    "mousePerturb4Tf" : "Xin Jin Lab Mouse KO Perturb-Seq, split per cell type",
    "comp" : "Xin Jin Lab Cell type composition change",
    "zebraAsdSizes" : "Hoffman Lab Brain Region Sizes",
    "perturbFishAstro" : "Farhi Lab Perturb-Fish Gene Expression Change in Astrocytes"
}
mousePerturbSql = 'select hgnc.symbol, m.gene, m.cellType, m.logFC, m.PValue, m.perturbation, m.target from mousePerturb4Tf as m, hgnc where hgnc.hgnc_id=m.hgnc_id order by perturbation;'
mousePerturbFields = [
    ("Perturbed Gene", 50),
    ("Mouse Effect Gene", 50),
    ("Cell Type", 50),
    ("LogFC", 50),
    ("pVal", 50),
    ("Guide", 50), 
    ("Perturbed Gene (mouse)", 50)
    ]
phenoNames = [x for x,y in mousePerturbFields]

compSql = "SELECT symbol, subcluster, PropRatio, FDR FROM comp, hgnc where comp.hgnc_id=hgnc.hgnc_id;" # "target" is the mouse gene
compFields = [
    ("Perturbed gene", 50),
    ("Cell Type", 50),
    ("Ratio change", 50),
    ("FDR", 50)
    ]

sizeSql = "select symbol, Forebrain, Optic_Tectum, Thalamus, Hypothalamus, Cerebellum, Hindbrain, Habenula, Posterior_Tuberculum,Mutant_Experiment_Sample from zebraAsdSizes as z, hgnc where hgnc.hgnc_id=z.hgnc_id;"
sizeFields = [
    ("KO gene (Human)", 50),
    ("Foreb Size", 50),
    ("Optic Tec Size", 50),
    ("Thamalus", 50),
    ("Hypoth", 50),
    ("Cereb", 50),
    ("Hindb", 50),
    ("Haben", 50),
    ("Poster Tuber", 50),
    ("Assay Name", 50),
    ]

perturbfishAstroSql = "SELECT symbol, gene, LFC, qVal from perturbFishAstro as p, hgnc where hgnc.hgnc_id=p.hgnc_id;"
perturbfishAstroFields = [
    ("Perturbed gene", 50),
    ("Changed Gene", 50),
    ("LFC", 50),
    ("qVal", 50),
    ]

allSql = "select * from clusterprop as c, perturbFishAstro as p, mousePerturb4Tf as m, zebraAsdSizes as z where c.hgnc_id = p.hgnc_id and p.hgnc_id = m.hgnc_id and m.hgnc_id = z.hgnc_id;"
allFields =  []
    

def debug(msg):
    if DEBUG:
        print(msg+"<br>")
        sys.stdout.flush()

def htmlHeader():
    " print start of page "
    print("""
<html>
<head>
<link rel="stylesheet" href="pure-min.css" integrity="sha384-X38yfunGUhNzHpBaEBsWLO+A0HDYOQi8ufWDkZ0k9e0eXz/tH3II7uKZ9msv++Ls" crossorigin="anonymous">

<style>

body {
    box-sizing: border-box; /* padding does not add to the size of the parent */
    font-family: helvetica;
    padding-left: 5px;
}

td { word-wrap: break-word; }

.pure-button { border-radius: 4px; }

.mainButton {
    font-size: 85%;
    padding: 3px 10px 3px 10px;
}

#mainContent {
    margin-left: 200px;
    position: relative;
    width: auto;
    overflow: auto;
    z-index: 1;
}

#facetBar {
    background-color: #EAEAEA;
    width: 200px;
    position: fixed;
    height: 100%;
    box-sizing: border-box;
    border: 1px solid #888;
}
#facetContent {
    padding: 3px;
}

#facetTitle {
    background-color: #CCC;
    color: black;
    font-weight: bold;
    padding-left: 5px;
}

#mainTable {
    margin-top: 10px;
    margin-left: 0px;
    table-layout: fixed;
}
#mainTable th {
    vertical-align: top;
    font-size: 12px;
    border-bottom: 1px solid #F0F0F0;
    padding: 2px;
}
#mainTable tr {
    vertical-align: top;
    border-bottom: 1px solid #F0F0F0;
    padding: 0;
}
#mainTable td {
    padding: 2px;
}
#tableSummary {
    margin-left: 8px;
}

.filterBox {
    border-top: 1px solid #888;
    margin-top: 3px;
    padding-top: 3px;
}

.numFiltBox {
    display: grid;
    grid-template-columns: 70px 130px;
}

.numFiltSpan {
    */float:right;*/
}
.numFiltInput {
    width: 70px;
}

.perPage {
    display: inline;
    margin-left: 20px;
    margin-right: 20px;

}

.sortIcon {
    font-size: 12px;
    margin-left:4px;
    cursor:pointer;
    color: #888;
}
.checkbox {
    margin-right: 5px;
}

#fieldList {
    display: grid;
    grid-template-columns: 40px 150px auto;
    width: 900px;
}

.fieldGroup {
    grid-column: span 3;
}
</style>

<title>SSPsygene Phenotypes</title>

</head>
<body>
<script>
</script>
""")

def htmlFooter():
    " print end of page "
    print("""
</body>
</html>
    """)

def printTableHeader(fieldWidths, fieldLabels):
    " print the <thead> part of the HTML table, size according to fieldWidths. Return dict fieldName -> index in row "

    #showIdx = [names.index(f[0]) for f in fieldWidths]
    #fieldToIdx = [names.index(f[0]) for f in fieldWidths]
    #assert(len(showIdx)==len(fieldWidths))

    print("<thead>")
    print("<tr>")
    sortField = cgiArgs.getvalue("sortby")
    for i, nameWidth in enumerate(fieldWidths):
        name, width = nameWidth
        label = fieldLabels.get(name, name)
        label = label.replace("_", " ")
        sortSym = "&#9660;"
        sortOrder = 'asc'
        addStyle = ""
        if sortField==name:
            sortSym = "&#9650"
            sortOrder = 'desc'
            addStyle = "color:white; background-color:#666;"
        print("<th style='width:%dpx;%s'>%s<span class='sortIcon' data-sortorder='%s' data-sortby='%s' >%s</span></th>" %
            (width, addStyle, label, sortOrder, name, sortSym))
    print("</tr>")
    print("</thead>")
    print("<tbody>")

    #return fieldToIdx

def valToHtml(field, val):
    if field=="field":
        val = val.replace("|", ", ")
    elif field=="hgncId":
        val = val
    elif field=="alias_name" or field=="alias_symbol":
        val = val.strip('"').replace("|", ", ")

    return val

def isLikeTrue(val):
    " return True is key is set in fieldStorage and one of 1 on or true "
    if val in ["1", "on", "true"]:
        return True
    return False

def hasNonAlpha(inStr):
    " return True is string contains non-alphanumeric characters and not an undersscore "
    for c in inStr:
        if not c.isalpha() and c!="_":
            return False
    return True

def parseRanges(rangesStr):
    ranges = []
    for rangeStr in rangesStr.split(","):
        parts = rangeStr.split(":")
        if len(parts)!=2:
            return None
        chrom, posRange = parts
        if hasNonAlpha(chrom):
            return None
        parts = posRange.split("-")
        if len(parts)!=2:
            return None
        start, end = parts
        start = int(start)
        end = int(end)
        ranges.append( (chrom, start, end) )
    return ranges

def cgiGetValsPrefix(prefix, valType, default=None):
    vals = []
    for key in cgiArgs.keys():
        if key.startswith(prefix):
            vals.append( (key, cgiGetVal(key, valType, default) ) )
    return vals

def cgiGetVal(name, valType, default=None):
    """ return value of cgi argument, make sure that input follows the specified format. 
    valType can be:
        - bool: 0,1, on,off or true or false,
        - ranges: chr1:0-1000,chr2:023-123
        - float
        - int
        - string: no brackets
    """
    val = cgiArgs.getvalue(name)
    if val is None:
        return default

    if valType=="bool":
        retVal = isLikeTrue(val)
    elif valType=="ranges":
        retVal = parseRanges(val)
        if retVal is None:
            print("Error: parameter %s has invalid format, must be comma-sep list of chrom:start-end, and chrom can only contain letters, numbers or underscore." % name)
    elif valType=="string":
        retVal = val.replace("<", "").replace(">", "").replace("'","").replace('"', '')
    elif valType=="float":
        retVal = float(val)
    elif valType=="int":
        retVal = int(val)
    else:
        assert(False)

    return retVal

def makeSql(args, joinTables, doCount=False, doPaging=False):
    " return full sql command to get main table "
    baseSql = "SELECT * FROM"
    if doCount:
        baseSql = "SELECT count(*) FROM"
    #sql = baseSql+" "+",".join(tables)
    sql = baseSql+" hgnc"

    if cgiGetVal("filter_tom", "bool") and not "tom" in showTables:
        showTables.append("tom")

    debug("Join on tables: %s" % joinTables)
    # make the JOIN ... parts
    joinList = []

    if "omim" in joinTables:
        joinList.append(("omim", "omim.omimId=hgnc.omim_id"))
    if "refseq_hg38" in joinTables:
        joinList.append(("refseq_hg38", "hgnc.refseq_accession=refseq_hg38.refseqId"))
    if "sfari" in joinTables:
        joinList.append(("sfari", "hgnc.hgnc_id=sfari.hgnc_id"))
    if "tom" in joinTables:
        joinList.append(("tom", "hgnc.hgnc_id=tom.hgnc_id"))
    if "singhFu" in joinTables:
        joinList.append(("singhFu", "hgnc.hgnc_id=singhFu.hgnc_id"))
    if "gnomad" in joinTables:
        joinList.append(("gnomad", "hgnc.hgnc_id=gnomad.hgnc_id"))

    joinStmts = []
    for table, fieldSpec in joinList:
        joinStmts.append("LEFT OUTER JOIN %s ON %s" % (table, fieldSpec))

    # make the WHERE parts
    whereList = []
    if cgiGetVal("filter_protein", "bool"):
        whereList.append('locus_group="protein-coding gene"')
    if cgiGetVal("filter_sfari", "bool"):
        whereList.append('gene_score is not NULL')
    if cgiGetVal("filter_syndromic", "bool"):
        whereList.append('syndromic<>"0"')

    for key, val in cgiGetValsPrefix("filt_", "float", None):
        if val!=None:
            fieldName = key.replace("filt_", "")
            op = cgiGetVal("op_"+fieldName, "string")
            if op=="gt":
                sqlOp = ">"
            elif op=="eq":
                sqlOp = "="
            elif op=="lt":
                sqlOp= "<"
            else:
                print("op_ for %s makes no sense" % fieldName)
                continue

            whereList.append(fieldName + " " + sqlOp + " " + str(val))


    filtRanges = cgiGetVal("filtRanges", "ranges")
    if filtRanges:
        rangeExprs = []
        for chrom, start, end in filtRanges:
            rangeExprs.append('chrom="%s" and start>=%d and end<=%d' % (chrom, start, end))
        combRangeExpr = " OR ".join(rangeExprs)
        whereList.append(combRangeExpr)

    # sorting
    orderList = []
    orderField = cgiArgs.getvalue("sortby")
    orderDir = cgiArgs.getvalue("sortorder", "asc")
    if orderField:
        # special hack to get the chrom pos sort order right
        if orderField=="chrom":
            orderField = "chrom,start,end"
        orderList.append("ORDER BY "+orderField)
        if orderDir!="asc":
            orderList.append("DESC")

    # limiting 
    limitStr = ""
    if doPaging:
        perPage = cgiGetVal("perPage", "int", 1000)
        page = cgiGetVal("page", "int", 1)
        fromLimit = (page-1) * perPage
        limitStr = " LIMIT %d OFFSET %d" % (perPage, fromLimit)

    # put everything together
    sql = sql+" "+" ".join(joinStmts)
    if len(whereList)!=0:
        sql = sql + " WHERE " + (" AND ".join(whereList))
    if len(orderList)!=0:
        sql = sql + " " + " ".join(orderList)
    sql = sql + limitStr

    debug("SQL: "+sql)
    return sql

def humReadable(num):
    if num > 1000000:
        return str(int(num/1000000))+"mbp"
    elif num > 1000:
        return str(int(num/1000))+"kbp"
    else:
        return num+"kbp"

def urlUcsc(db, chrom, start, end):
    posStr = "%s:%s-%s" % (chrom, start, end)
    return "https://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s" % (db, posStr)

def notEmpty(x):
    " sqlite returns None for empty fields, tsv files use empty strings "
    return x is not None and x!=""

def rowToHtml(fields, row):
    htmlRow = []

    for name, val in zip(fields, row):
        htmlVal = valToHtml(name, val)

        if name=="symbol":
            longName = row[fieldToIdx["name"]]
            htmlVal = "%s<br><small>%s</small>" % (htmlVal, longName)
        elif name=="chrom":
            if notEmpty(val):
                chrom = val
                start = row[fieldToIdx["start"]]
                end = row[fieldToIdx["end"]]
                hg38Url = urlUcsc("hg38", chrom, start, end)
                htmlVal = "<small><a target=_blank href='%s'>%s:%s</a></small><br>" % \
                    (hg38Url, chrom, humReadable(int(start)))
        elif name=="pheno":
            if val is not None:
                htmlVal = val.replace("|", "; ");
        if val is None:
            htmlVal = ""

        if isinstance(val, float):
            if val!=0.0:
                htmlVal = "%.2f" % val
            else:
                htmlVal = "0"

        htmlRow.append( "<td>%s</td>" % htmlVal )
    return "\n".join(htmlRow)

def printTiming(time):
    timeDiff = timeit.default_timer() - time
    if measureTiming:
        print("Time: %s" % timeDiff)

def printCheckbox(name, label, addBr=False, isChecked=False):
    checkStr = ""
    if cgiGetVal(name, "bool") or isChecked:
        checkStr = " checked"
    print("<input name='%s' class='checkbox' type='checkbox'%s>%s</input>" % (name, checkStr, label))
    if addBr:
        print("<br>")

def printNumberFilter(name, label):
    defValue = cgiGetVal("filt_"+name, "float", "")
    defStr = str(defValue)
    print("<div>")
    print(label)
    print("</div>")

    print('<div class="numFiltSpan">')

    opVal = cgiGetVal("op_"+name, "string", None)
    defVals = {}
    if opVal is not None:
        defVals[opVal] = 'selected="selected"'

    print('<select name="op_%s">' % name)
    print('  <option value="lt" %s>&lt;</option>' % defVals.get("lt", ""))
    print('  <option value="eq" %s>=</option>' % defVals.get("eq", ""))
    print('  <option value="gt" %s>&gt;</option>' % defVals.get("gt", ""))
    print('</select>')

    print('<input class="numFiltInput" step="any" name="filt_%s" type="number" value="%s" />' % (name, defStr))
    print("</div>")
    #print("</div>") # numFiltBox

def getCursor():
    # get count of data rows
    conn = sqlite3.connect("pheno.db")
    cur = conn.cursor()
    return cur

def getRowCount(cur):
    # get data rows
    time = timeit.default_timer()
    sql = makeSql(cgiArgs, showTables, doCount=True)
    cur.execute(sql)
    printTiming(time)
    count = cur.fetchone()[0]
    return count

def getRowsTsv():
    " return the tsv rows "
    ifh = open("genes.tsv")
    headers = None
    rows = []
    for line in ifh:
        row = line.rstrip("\n").split("\t")
        if headers is None:
            headers = row
        else:
            rows.append(row)
    return rows, headers

def getRowsSqlite(cur):
    " return all rows of the current big table "
    time = timeit.default_timer()

    sql = makeSql(cgiArgs, showTables, doPaging=True)
    cur.execute(sql)

    colNames = [description[0] for description in cur.description]
    if DEBUG:
        print("Sql fields: %s" % colNames)
    printTiming(time)

    return cur.fetchall(), colNames

def downloadTable():
    #for row in parseTsv("hgnc_complete_set.txt"):
    cur = getCursor()
    rows, colNames = getRowsSqlite(cur)

    print("Content-type: text/tab-separated-values")
    print('Content-Disposition: attachment; filename="geneTable.csv"')
    print()

    print("\t".join(colNames))
    for row in rows:
        row = [str(x) for x in row]
        print("\t".join(row))

def parseFieldInfo():
    " parse fields.tsv "
    fields = []
    for line in open("fields.tsv"):
        if line.startswith ("#"):
            headers = line.lstrip("#").rstrip("\n").split("\t")
            continue

        row = line.rstrip("\n").split("\t")
        name, group, fieldType = row[:3]
        shortLabel = ""
        longLabel = ""
        if len(row)>3:
            shortLabel = row[3]
            longLabel = shortLabel
        if len(row)>4:
            longLabel = row[4]
        fields.append( (name, group, fieldType, shortLabel, longLabel) )
    return fields

def parseFieldInfo():
    " parse fields.tsv "
    fields = []
    for line in open("fields.tsv"):
        if line.startswith ("#"):
            headers = line.lstrip("#").rstrip("\n").split("\t")
            continue

        row = line.rstrip("\n").split("\t")
        name, group, fieldType = row[:3]
        desc = ""
        if len(row)>3:
            desc = row[3]
        fields.append( (name, group, fieldType, desc) )
    return fields

def pickFieldsPage():
    ""
    print("<h3>Select fields shown on table:</h3>")
    fields, groupLabels = parseFieldInfo()

    fieldsByGroup = defaultdict(list)
    groups = []
    for name, group, fieldType, desc in fields:
        fieldsByGroup[group].append( (name, fieldType, desc) )
        if group not in fieldsByGroup:
            groups.append(group)

    currFieldStr = cgiGetVal("fields", "string")
    if currFieldStr:
        currFields = currFieldStr.split(",")
    else:
        currFields = []

    for group in groups:
        shortLabel, longLabel = groupLabels[group]
        print(f"<a href='#{group}'>{shortLabel} - {longLabel}</a><br>")

    print("<form>")

    print('<input id="submitForm" type="submit" value="Save selection and return to table">')

    print("<div id='fieldList'>")

    for group, fieldInfos in fieldsByGroup.items():
        print(f"<h4 id='{group}' class='fieldGroup'>{group}</h4>")
        for fieldInfo in fieldInfos:
            name, fieldType, desc = fieldInfo
            label = "<div>"+name+"</div><div>"+desc+"</div>"
            printCheckbox(name, label, isChecked=(name in currFields))
    print("</div>") # fieldList

    print("<p>")
    print('<input id="submitForm" type="submit" value="Save selection and return to table">')
    print("</form>")

def printSelect(name, options, selVal):
    print('<select name="%s">' % name)
    for val, label in options:
        if val==selVal:
            selStr = " selected='selected'"
        else:
            selStr = ""
        print('  <option value="%s" %s>%s</option>' % (val, selStr, label))
    print('</select>')

def filterRowsTsv(rows, colNames):
    " "
    newRows = []
    for row in rows:
        newRows.append(row)
        
    return newRows

def getFields(res):
    " return list of fields given sql result "
    desc = res.description
    fields = [x[0] for x in desc]
    return fields
    
def tableToDict(cur, geneRowDict, fieldDict, table):
    " get all rows from table table and add to geneRowDict on first field "
    res = cur.execute("SELECT h.symbol, t.* from %s as t, hgnc as h where h.hgnc_id=t.hgnc_id" % table)
    fieldNames = getFields(res)
    fieldDict[table] = fieldNames[1:]

    for row in res.fetchall():
        sym = row[0]
        if sym not in geneRowDict:
            geneRowDict[sym] = {}
        if table not in geneRowDict[sym]:
            geneRowDict[sym][table] = []
        geneRowDict[sym][table].append( row[1:] )

    return geneRowDict, fieldDict

def printMerge(geneRowDict, fieldDict):
    " output the big merged gene symbol dict "
    print("<h2>All phenotypes for all perturbed/knockout genes</h2>")

    for sym, tableRows in geneRowDict.items():
        tableNames = list(tableRows.keys())
        tableLabels = [tableShort[t] for t in tableNames]

        print("<a href='#%s'>%s</a> - %d assays - %s<br>" % (sym, sym, len(tableRows), " - ".join(tableLabels)))
    print("<hr>")

    for sym, tableRows in geneRowDict.items():
        print(f"<h3 id={sym}>{sym}</h3>")
        print("<ul>")
        for table, rows in tableRows.items():
            tableDesc = tableDescs[table]
            print(f"<h4>{tableDesc}</h4>")
            print("<table class='pure-table'>")
            
            print("<thead><tr>")
            fieldNames = fieldDict[table]
            for fieldName in fieldNames:
                if fieldName=="hgnc_id":
                    continue
                print(f"<th>{fieldName}</th>")
            print("</tr></thead>")

            for row in rows:
                print("<tr>")
                for fieldName, val in zip(fieldNames, row):
                    if fieldName=="hgnc_id":
                        continue
                    print(f"<td>{val}</td>")
                print("</tr>")
            print("</table>")
        print("</ul>")

def htmlMiddle(args):
    " print html middle part "

    print("<h1>SSPsyGene Phenotype Knowledgebase</h1>")
    # https://docs.google.com/spreadsheets/d/1q5pjIyDGnz4L3DNh7npbKE2TCL_IGgJU0njEEEY2n-U/edit#gid=0

    print("<div id='mainContent'>")

    page = args.getvalue("page")

    if page is None:
        print("<h4>Individual assays:</h4>")
        print('<p><a href="?page=deg">Xin Jin: Mouse Perturb-Seq Gene/Gene expression changes</a></p>')
        print('<p><a href="?page=comp">Xin Jin: Mouse Perturb-Seq Gene/Cell type composition changes</a></p>')
        print('<p><a href="?page=sizes">Ellen Hoffman: Zebrafish Brain Region Sizes - Gene/Brain Size </a></p>')
        print('<p><a href="?page=perturbFishAstr">Sami Farhi: Human Astrocyte Perturb-Fish Expr Changes - Gene/Gene </a></p>')
        print("<h4>All assays:</h4>")
        print('<p><a href="?page=all">Integrated assays</a></p>')
    else:
        cur = getCursor()
        print("<p><a href='show.py'>Back</a></p>")

        if page=="deg":
            sql = mousePerturbSql
            title = "Gene expression changes"
            fields = mousePerturbFields

        elif page=="comp":
            sql = compSql
            fields = compFields
            title = "Cell type composition changes"

        elif page=="sizes":
            sql = sizeSql
            fields = sizeFields
            title = "Zebrafish brain region sizes"

        elif page=="perturbFishAstr":
            sql = perturbfishAstroSql
            fields = perturbfishAstroFields
            title = "Perturb-fish expression changes"

        elif page=="all":
            tables = ["mousePerturb4Tf", "comp", "zebraAsdSizes", "perturbFishAstro"]
            geneRowDict = {} # sym -> table -> list of rows
            fieldDict = {} # table -> list of fieldnames
            for t in tables:
                geneRowDict, fieldDict = tableToDict(cur, geneRowDict, fieldDict, t)

            printMerge(geneRowDict, fieldDict)
            sys.exit(0)

        else:
            print("invalid value for 'page'")

        print("<table style='' id='mainTable' class='pure-table'>")
        print("<h4>%s</h4>" % title)
        printTableHeader(fields, {})

        res = cur.execute(sql)
        
        for row in res.fetchall():
            print("<tr>")
            print( rowToHtml(fields, row) )
            print("</tr>")

        print("</tbody></table>")

    print("</div>") # mainContent

def getCookies(cookieDb):
    """ return db cookie value. Called usually before headers are printed.
    Default value can be sent as argument.
    """
    cookie_string = os.environ.get('HTTP_COOKIE')
    if cookie_string:
        cookie = Cookie.SimpleCookie()
        cookie.load(cookie_string)
        if "db" in cookie:
            cookieDb = cookie["db"].value
    return cookieDb

def setCookies(db):
    """
    Send a cookie header to set the "db" cookie to the value specified
    """
    cookie = Cookie.SimpleCookie()
    cookie['db'] = db
    cookie['db']['expires'] = 30 * 24 * 60 * 60
    cookie['db']['comment'] = 'holds the last hgMirror database'
    print(cookie)

def main():
    #defaultDb = getCookies(defaultDb)
    args = cgi.FieldStorage()
    global cgiArgs
    cgiArgs = args

    if "download" in args:
        downloadTable()
        sys.exit(0)

    print("Content-type: text/html\n")

    #if "measureTiming" in args:
        #global measureTiming
        #measureTiming = True
    if "debug" in args:
        global DEBUG
        DEBUG=True

    htmlHeader()
    htmlMiddle(args)
    htmlFooter()

main()

