import json,sys,re
nb=json.load(open(sys.argv[1])); kw=re.compile(sys.argv[2],re.I); maxc=int(sys.argv[3])
for i,c in enumerate(nb['cells']):
    src=''.join(c.get('source',[]))
    if c['cell_type']=='markdown':
        if kw.search(src): print(f'--- md {i}:',src[:600].replace('\n',' | '))
        continue
    outs=[]
    for o in c.get('outputs',[]):
        t=''.join(o.get('text',[])) if 'text' in o else ''.join(o.get('data',{}).get('text/plain',[]))
        if t and not t.startswith('<Figure'): outs.append(t)
    txt='\n'.join(outs)
    if kw.search(src) or kw.search(txt):
        print(f'=== code {i}:',src[:500].replace('\n',' | ')); 
        if txt: print('  OUT:',txt[:maxc])
