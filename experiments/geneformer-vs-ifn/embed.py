"""Frozen Geneformer-V2-104M cell embeddings. Usage: python -I embed.py <gf_dir> <q3_dir>"""
import sys, os, pickle, numpy as np, torch
from transformers import BertModel
gfd, d = sys.argv[1:3]
torch.set_num_threads(2)
m = BertModel.from_pretrained(f'{gfd}/Geneformer-V2-104M', attn_implementation='sdpa').eval()
seqs = pickle.load(open(f'{d}/tokens.pkl', 'rb'))
out = f'{d}/emb_partial.npy'
E = np.load(out) if os.path.exists(out) else np.full((len(seqs), 768), np.nan, dtype=np.float32)
order = np.argsort([len(s) for s in seqs]); B = 4
for k in range(0, len(order), B):
    b = order[k:k + B]
    if not np.isnan(E[b]).any(): continue
    Lm = max(len(seqs[i]) for i in b)
    ids = torch.zeros((len(b), Lm), dtype=torch.long); mask = torch.zeros_like(ids)
    for r, i in enumerate(b):
        ids[r, :len(seqs[i])] = torch.tensor(seqs[i]); mask[r, :len(seqs[i])] = 1
    with torch.no_grad():
        h = m(input_ids=ids, attention_mask=mask, output_hidden_states=True).hidden_states[-2]
    for r, i in enumerate(b):
        n = len(seqs[i]); E[i] = h[r, 1:n - 1].mean(0).numpy()   # exclude <cls> and <eos>
    if (k // B) % 25 == 0:
        np.save(out, E); print(k + len(b), 'done', flush=True)
np.save(out, E); np.save(f'{d}/emb.npy', E); print('finished', flush=True)
