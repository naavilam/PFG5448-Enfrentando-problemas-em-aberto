# Ink Boat — Marangoni Propulsion

Aplicação Streamlit centrada na **força e velocidade de propulsão** causadas por gradientes de tensão superficial.

O núcleo numérico é o mesmo da versão `inkboat_coupled`: transporte superficial efetivo + solução espectral de Stokes para uma subfase profunda com superfície plana + força capilar na linha de contato + arrasto aproximado. **Não é um solver trifásico completo.** As tensões tinta–ar e tinta–água são exibidas no coeficiente S, mas não entram diretamente na força propulsiva. Válido apenas sob as hipóteses declaradas, especialmente Re pequeno.

## Executar

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Publicar

No Streamlit Community Cloud, escolher o caminho da subpasta contendo `app.py`. Atualizações na branch publicada são aplicadas após `git push`.
