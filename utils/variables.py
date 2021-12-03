#                combination name : [combination's gene list:list, image size:tuple]
combinations = {'immunomodulation_markers': [['CCL5', 'CCR5', 'CCR3'], (10, 6)],
                'differentiation_markers': [['GFAP', 'BETATUBULIN'], (10, 6)],
                'cannaboid_receptors_markers': [['GPR55', 'CNR1', 'TRPV1'], (10, 6)],
                'gbm_stem_cells_markers': [['CD9', 'SOX2', 'CD44', 'CD15', 'ID1', 'OCT4',
                                            'TRIM28', 'TUFM', 'OLIG2', 'NOTCH', 'ALYREF'], (15, 11)],
                'emt_markers': [['CDH1', 'VIM', 'SNAI1', 'NF-KB', 'CD44', 'STAT3', 'CHI3L1', 'NOTCH'], (15, 11)],
                'mix_tcga_subtype_markers': [['ACSBG1', 'KCNF1', 'S100A4', 'NF-KB', 'P2RX7', 'STMN4', 'SOX10', 'ERBB3',
                                              'OLIG2', 'NOTCH', 'COL1A2', 'COL1A', 'TGFB1', 'THBS1', 'DAB2'], (20, 16)],
                'classical_tcga_subtype_markers': [['ACSBG1', 'KCNF1', 'S100A4', 'NF-KB'], (10, 6)],
                'proneural_tcga_subtype_markers': [['P2RX7', 'STMN4', 'SOX10', 'ERBB3', 'OLIG2', 'NOTCH'], (10, 6)],
                'mesenchymal_tcga_subtype_markers': [['COL1A2', 'COL1A', 'TGFB1', 'THBS1', 'DAB2'],
                                                     (13, 9)]}

subtypes = {'classical': ['ACSBG1', 'KCNF1', 'S100A4', 'NF-KB'],
            'proneural': ['P2RX7', 'STMN4', 'SOX10', 'ERBB3', 'OLIG2', 'NOTCH'],
            'mesenchymal': ['COL1A2', 'COL1A', 'TGFB1', 'THBS1', 'DAB2']}

subtype_km = {'CL': ['ACSBG1', 'KCNF1', 'S100A4', 'NF-KB'],
              'PN': ['P2RX7', 'STMN4', 'SOX10', 'ERBB3', 'OLIG2', 'NOTCH'],
              'MES': ['COL1A2', 'COL1A', 'TGFB1', 'THBS1', 'DAB2'],
              'MIX': ['COL1A2', 'COL1A', 'TGFB1', 'THBS1', 'DAB2', 'S100A4', 'P2RX7', 'STMN4', 'SOX10',
                      'ERBB3', 'ACSBG1', 'OLIG2', 'NOTCH', 'NF-KB', 'KCNF1']}

mix_tcga = ['ACSBG1', 'KCNF1', 'S100A4', 'NF-KB', 'P2RX7', 'STMN4', 'SOX10', 'ERBB3',
            'OLIG2', 'NOTCH', 'COL1A2', 'COL1A', 'TGFB1', 'THBS1', 'DAB2']

comb_groups = ['immunomodulation_markers', 'differentiation_markers', 'cannaboid_receptors_markers',
               'gbm_stem_cells_markers',
               'emt_markers', 'mix_tcga_subtype_markers', 'classical_tcga_subtype_markers',
               'proneural_tcga_subtype_markers',
               'mesenchymal_tcga_subtype_markers']

groups = [('immunomodulation_markers', 'differentiation_markers', (10, 6)),
          ('immunomodulation_markers', 'cannaboid_receptors_markers', (10, 6)),
          ('immunomodulation_markers', 'gbm_stem_cells_markers', (15, 11)),
          ('immunomodulation_markers', 'emt_markers', (12, 8)),
          ('immunomodulation_markers', 'mix_tcga_subtype_markers', (16, 12)),
          ('immunomodulation_markers', 'classical_tcga_subtype_markers', (10, 6)),
          ('immunomodulation_markers', 'proneural_tcga_subtype_markers', (10, 6)),
          ('immunomodulation_markers', 'mesenchymal_tcga_subtype_markers', (10, 6)),
          ('differentiation_markers', 'cannaboid_receptors_markers', (10, 6)),
          ('differentiation_markers', 'gbm_stem_cells_markers', (15, 11)),
          ('differentiation_markers', 'emt_markers', (11, 7)),
          ('differentiation_markers', 'mix_tcga_subtype_markers', (17, 13)),
          ('differentiation_markers', 'classical_tcga_subtype_markers', (10, 6)),
          ('differentiation_markers', 'proneural_tcga_subtype_markers', (10, 6)),
          ('differentiation_markers', 'mesenchymal_tcga_subtype_markers', (10, 6)),
          ('cannaboid_receptors_markers', 'gbm_stem_cells_markers', (15, 11)),
          ('cannaboid_receptors_markers', 'emt_markers', (13, 9)),
          ('cannaboid_receptors_markers', 'mix_tcga_subtype_markers', (17, 13)),
          ('cannaboid_receptors_markers', 'classical_tcga_subtype_markers', (10, 6)),
          ('cannaboid_receptors_markers', 'proneural_tcga_subtype_markers', (10, 6)),
          ('cannaboid_receptors_markers', 'mesenchymal_tcga_subtype_markers', (11, 7)),
          ('gbm_stem_cells_markers', 'emt_markers', (18, 14)),
          ('gbm_stem_cells_markers', 'mix_tcga_subtype_markers', (22, 18)),
          ('gbm_stem_cells_markers', 'classical_tcga_subtype_markers', (15, 11)),
          ('gbm_stem_cells_markers', 'proneural_tcga_subtype_markers', (16, 12)),
          ('gbm_stem_cells_markers', 'mesenchymal_tcga_subtype_markers', (18, 14)),
          ('emt_markers', 'mix_tcga_subtype_markers', (18, 14)),
          ('emt_markers', 'classical_tcga_subtype_markers', (12, 8)),
          ('emt_markers', 'proneural_tcga_subtype_markers', (12, 8)),
          ('emt_markers', 'mesenchymal_tcga_subtype_markers', (15, 11)),
          ('classical_tcga_subtype_markers', 'proneural_tcga_subtype_markers', (12, 8)),
          ('classical_tcga_subtype_markers', 'mesenchymal_tcga_subtype_markers', (10, 6)),
          ('proneural_tcga_subtype_markers', 'mesenchymal_tcga_subtype_markers', (10, 6))]

subtype_target_names = {0: 'CL', 1: 'MES', 2: 'MIX', 3: 'PN'}

genes = ['TUFM', 'TRPV1', 'TRIM28', 'THBS1', 'TGFB1', 'STMN4', 'SOX10', 'S100A4', 'P2RX7', 'KCNF1', 'GPR55', 'ERBB3',
         'DAB2', 'CST7', 'COL1A2', 'COL1A', 'CNR1', 'CD9', 'CCR5', 'CCL5', 'ACSBG1', 'VIM', 'STAT3', 'SOX2', 'SNAI1',
         'PROM1', 'OLIG2', 'OCT4', 'NOTCH', 'NF-KB', 'ID1', 'GFAP', 'CHI3L1', 'CEBPA', 'CDH1', 'CD44', 'CD15', 'CCR3',
         'BETATUBULIN', 'ALYREF']

biomarkers_for_subtype_definition = ['COL1A2', 'COL1A', 'TGFB1', 'THBS1', 'DAB2', 'P2RX7', 'STMN4', 'SOX10',
                                     'ERBB3', 'ACSBG1', 'OLIG2', 'NOTCH', 'NF-KB', 'S100A4', 'KCNF1']

