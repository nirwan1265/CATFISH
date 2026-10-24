import re, os

MAG = os.environ['HOME'] + '/mnt/MAGCAT'
THE = os.environ['HOME'] + '/mnt/Thesis_NirwanTandukar_GeneticsGenomics/ncsuthesis-0.6'

oup = open(MAG + '/oup-authoring-template.tex', encoding='utf-8').read()
ch  = open(THE + '/Chapter-5/Chapter-5.tex', encoding='utf-8').read()

def cut(txt, start_pat, end_pat, inc_start=True, inc_end=False):
    s = re.search(start_pat, txt); e = re.search(end_pat, txt[s.end():])
    a = s.start() if inc_start else s.end()
    b = s.end() + (e.end() if inc_end else e.start())
    return txt[a:b]

# ---------- 1. front matter ----------
front = oup[:re.search(r'\\section\{1\.1 Introduction\}', oup).start()]
# add missing table packages + fix graphicspath
front = front.replace(r'\usepackage{algorithm}',
    '% give LaTeX room to place the full-width floats\n\\renewcommand{\\topfraction}{0.92}\n\\renewcommand{\\bottomfraction}{0.7}\n\\renewcommand{\\textfraction}{0.08}\n\\renewcommand{\\floatpagefraction}{0.75}\n\\renewcommand{\\dbltopfraction}{0.92}\n\\renewcommand{\\dblfloatpagefraction}{0.75}\n\\setcounter{topnumber}{3}\n\\setcounter{totalnumber}{5}\n\\setcounter{dbltopnumber}{3}\n\\usepackage{booktabs}\n\\usepackage{makecell}\n\\usepackage{multirow}\n\\usepackage{array}\n\\usepackage{tabularx}\n\\usepackage{pdflscape}\n\\usepackage{algorithm}')
front = front.replace(r'\graphicspath{{Fig/}}', '%\\graphicspath{{Fig/}}  % paths are given explicitly below')

# ---------- 2. body from the thesis chapter ----------
body = cut(ch, r'\\section\{Introduction\}', r'\\section\{Data availability\}')

# strip thesis-only scaffolding
body = re.sub(r'\\label\{ch5-[^}]*\}', '', body)
body = body.replace('Figures/Chapter5/main/', 'Final/final_main_fig/')
body = body.replace('Figures/Chapter5/supplementary/', 'Final/final_sup/supp_figs/')

# section/subsection renumbering  (longest/most specific first)
ren = [
 (r'\section{Introduction}', r'\section{1 Introduction}'),
 (r'\subsection{Pathway signal archetypes}', r'\subsection*{1.1 Pathway signal archetypes}'),
 (r'\section{Methods}', r'\section{2 Methods}'),
 (r'\subsection{Overview of CATFISH and its component tests}', r'\subsection*{2.1 Overview of CATFISH and its component tests}'),
 (r'\subsection{Dependence structure of the tests}', r'\subsection*{2.2 Dependence structure of the tests}'),
 (r'\subsection{LD-aware MVN calibration}', r'\subsection*{2.3 LD-aware MVN calibration}'),
 (r'\subsubsection{Strategy 1: calibrate components, then combine.}', r'\subsubsection*{2.3.1 Strategy 1: calibrate components, then combine}'),
 (r'\subsubsection{Strategy 2: combine, then calibrate.}', r'\subsubsection*{2.3.2 Strategy 2: combine, then calibrate}'),
 (r'\subsection{Candidate-gene prioritization}', r'\subsection*{2.4 Candidate-gene prioritization}'),
 (r'\subsection{Simulation-based evaluation of CATFISH}', r'\subsection*{2.5 Simulation-based evaluation of CATFISH}'),
 (r'\section{Results}', r'\section{3 Results}'),
 (r'\subsection{Archetype-specific response profiles}', r'\subsection*{3.1 Archetype-specific response profiles}'),
 (r'\subsection{LD-aware null calibration and sensitivity to incomplete correlations}', r'\subsection*{3.2 LD-aware null calibration and sensitivity to incomplete correlations}'),
 (r'\subsection{Stress tests and diagnostic perturbations}', r'\subsection*{3.3 Stress tests and diagnostic perturbations}'),
 (r'\subsection{Integrated GWAS, MAGMA, and CATFISH analyses identify pathways associated with stem volume in sorghum}', r'\subsection*{3.4 Integrated GWAS, MAGMA, and CATFISH analyses identify pathways associated with stem volume in sorghum}'),
 (r'\subsection{CATFISH component tests reveal complementary pathway signals and refine candidate-gene prioritization for sorghum stem volume}', r'\subsection*{3.5 CATFISH component tests reveal complementary pathway signals and refine candidate-gene prioritization for sorghum stem volume}'),
 (r'\subsection{Phenotype-permutation null calibration of CATFISH in SAP stem volume}', r'\subsection*{3.6 Phenotype-permutation null calibration of CATFISH in SAP stem volume}'),
 (r'\section{Discussion}', r'\section{4 Discussion}'),
 (r'\subsection{Interpretation of archetype-specific response profiles}', r'\subsection*{4.1 Interpretation of archetype-specific response profiles}'),
 (r'\subsection{Dependence-aware calibration is central to CATFISH robustness}', r'\subsection*{4.2 Dependence-aware calibration is central to CATFISH robustness}'),
 (r'\subsection{4.3 Pathway-aware prioritization refines candidate-gene interpretation for sorghum stem volume}', r'\subsection*{4.3 Pathway-aware prioritization refines candidate-gene interpretation for sorghum stem volume}'),
 (r'\subsection{4.4 Empirical phenotype-permutation calibration highlights both the value and limits of MVN correction}', r'\subsection*{4.4 Empirical phenotype-permutation calibration highlights both the value and limits of MVN correction}'),
 (r'\section{Conclusion}', r'\section{5 Conclusion}'),
]
miss = []
for a, b in ren:
    if a in body: body = body.replace(a, b, 1)
    else: miss.append(a)

# ---- float conversion: thesis is one-column, OUP class is two-column ----
# full-width floats must be starred (figure*/table*); landscape pages are
# replaced by full-width starred floats, which is the OUP idiom.
body = re.sub(r'\\(begin|end)\{landscape\}\s*\n', '', body)

def star_floats(txt):
    lines = txt.split('\n')
    out, stack = [], []
    for ln in lines:
        mb = re.match(r'(\s*)\\begin\{(figure|table)\}(\[[^\]]*\])?', ln)
        me = re.match(r'(\s*)\\end\{(figure|table)\}', ln)
        if mb:
            ind, env, opt = mb.group(1), mb.group(2), mb.group(3) or ''
            place = '[!tp]' if env == 'table' else '[!t]'
            out.append(f'{ind}\\begin{{{env}*}}{place}' + ln[mb.end():])
            stack.append(env)
        elif me and stack and stack[-1] == me.group(2):
            env = stack.pop()
            out.append(f'{me.group(1)}\\end{{{env}*}}' + ln[me.end():])
        else:
            out.append(ln)
    return '\n'.join(out)

body = star_floats(body)

# appendix cross-refs
body = body.replace('Appendix 3.7', 'Appendix~A3.7').replace('Appendix A.5', 'Appendix~A5')
body = body.replace('is in the Appendix.', 'is in Appendix~A2.')

# ---------- 3. back matter: sections 6-10 ----------
back = cut(oup, r'\\section\{6 Conflicts of interest\}', r'\\begin\{thebibliography\}')
# swap in the chapter's real data-availability text
back = re.sub(r'(\\section\{8 Data availability\}\n).*?(\n\\section\{9)',
  r'\1\nThe R package and all code used to generate the figures and run the simulations underlying '
  r'this article are available on GitHub at \\href{https://github.com/nirwan1265/CATFISH}{github.com/nirwan1265/CATFISH}.\n\2',
  back, flags=re.S)

# ---------- 4. appendices, up to the Supplementary Figures section ----------
app = cut(oup, r'\\begin\{appendices\}', r'\\section\{Supplementary Figures\}')

# fix a stale citation key carried over from the old draft
app = app.replace('deLeeuw2015MAGMA', 'deLeeuw2015')

# ---------- 5. supplementary figures from the chapter ----------
supp = cut(ch, r'\\section\{Supplementary figures\}', r'\\section\{Supplementary tables\}')
supp = supp.replace('Figures/Chapter5/supplementary/', 'Final/final_sup/supp_figs/')
supp = re.sub(r'\\renewcommand\{\\thefigure\}[^\n]*\n', '', supp)
supp = re.sub(r'\\setcounter\{figure\}\{0\}', '', supp)
supp = re.sub(r'%%+ ETD:[^\n]*\n(%[^\n]*\n)*', '', supp)
supp = supp.replace(r'\section{Supplementary figures}', r'\section{Supplementary Figures}')
# number these S1..S4 so they read as supplementary, not as Figures 11-14
supp = supp.replace(r'\section{Supplementary Figures}',
    '\\section{Supplementary Figures}\n\\setcounter{figure}{0}\n\\renewcommand{\\thefigure}{S\\arabic{figure}}')

out = '\n'.join([
  front,
  '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
  '%% MAIN TEXT -- ported from thesis Chapter-5 (updated 2026-10-23)',
  '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
  '', body, '',
  back,
  '\\bibliographystyle{oup-abbrvnat}',
  '\\bibliography{Chapter5-CATFISH-refs}',
  '', '\\clearpage', '\\newpage', '',
  app,
  supp,
  '\\end{appendices}',
  '\\end{document}',
])

open(MAG + '/CATFISH_manuscript.tex', 'w', encoding='utf-8').write(out)
print("wrote CATFISH_manuscript.tex  lines:", out.count('\n'), " chars:", len(out))
print("unmatched section patterns:", miss if miss else "none")
