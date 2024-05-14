import numpy as np


antenna_coords=np.array([[-0.0271,-3.6527,6370342.6829],
[-0.0206,-3.5359,6370341.2279],
[-0.0332,-3.6390,6370342.7749],
[0.0203,-3.6260,6370342.2497],
[-0.0403,-3.5588,6370341.3543],
[0.0476,-3.6042,6370341.2712],
[-0.0784,-3.6016,6370342.0960],
[-0.1011,-3.6007,6370342.0313],
[-0.0102,-3.5804,6370341.2182],
[0.0450,-3.5876,6370341.5134],
[-0.0074,-3.6157,6370342.1803],
[-0.0276,-3.6597,6370342.4935],
[-0.0243,-3.6133,6370379.6591],
[-0.0530,-3.6088,6370341.6511],
[0.0191,-3.5803,6370341.5564],
[-0.0240,-3.5930,6370379.3082],
[-0.0598,-3.6448,6370380.2054],
[0.0111,-3.5959,6370379.3583],
[0.0694,-3.6331,6370341.4418],
[-0.0740,-3.6221,6370379.8124],
[-0.0756,-3.5738,6370342.0741],
[0.0438,-3.6361,6370342.1351],
[0.0471,-3.5622,6370341.8037],
[-0.0577,-3.5470,6370341.5590],
[0.0063,-3.6389,6370380.1026],
[0.0082,-3.5605,6370378.7444],
[-0.0459,-3.5886,6370379.2304],
[-0.0155,-3.5623,6370378.7748],
[-0.0900,-3.6331,6370380.0033],
[-0.1830,-3.8470,6370339.4602],
[0.0932,-3.8200,6370341.9617],
[-0.1955,-3.5079,6370340.2910],
[-0.1065,-3.7563,6370339.8583],
[-0.0273,-3.7071,6370343.3455],
[-0.0875,-3.6574,6370341.6531],
[-0.1725,-3.6449,6370337.7768],
[-0.0379,-3.7179,6370342.2231],
[-0.0308,-3.7827,6370342.6859],
[-0.1357,-3.6311,6370340.2580],
[-0.1037,-3.7201,6370340.5113],
[-0.1249,-3.7460,6370339.5499],
[-0.1229,-3.7029,6370339.1827],
[-0.1209,-3.6374,6370339.7073],
[-0.0571,-3.6855,6370342.6008],
[-0.1481,-3.6881,6370338.0962],
[-0.1683,-3.6624,6370380.5109],
[-0.1108,-3.6646,6370380.5495],
[-0.0802,-3.6839,6370380.8827],
[-0.0728,-3.7837,6370341.1838],
[0.0865,-3.5052,6370341.4351],
[-0.0559,-3.7048,6370342.4154],
[-0.0536,-3.7589,6370343.0526],
[-0.0515,-3.7395,6370342.8762],
[-0.0781,-3.7382,6370340.7449],
[-0.1542,-3.6247,6370379.8564],
[-0.0286,-3.7538,6370342.7145],
[-0.1152,-3.6777,6370380.7760],
[-0.1484,-3.7252,6370338.3885],
[-0.0594,-3.6635,6370380.5288],
[-0.0441,-3.6869,6370380.9346],
[-0.1350,-3.6563,6370380.4049],
[-0.0843,-3.7057,6370381.2613],
[0.0795,-3.8523,6370342.6922],
[-0.1453,-3.5157,6370341.7479],
[0.0339,-3.6960,6370341.8230],
[0.0345,-3.6942,6370342.0320],
[0.0337,-3.6994,6370341.9521],
[0.0358,-3.6967,6370342.3845],
[0.0352,-3.7013,6370342.0656],
[0.0376,-3.7003,6370342.3067],
[0.0350,-3.7037,6370341.8566],
[0.0390,-3.7015,6370342.3180],
[0.0371,-3.7046,6370341.9318],
[0.0334,-3.7041,6370341.8641],
[0.0369,-3.7058,6370341.9029],
[0.0353,-3.7024,6370341.9629],
[0.0308,-3.6987,6370341.8089],
[0.0323,-3.6976,6370341.5400],
[0.0317,-3.7010,6370341.8590],
[0.0284,-3.7006,6370342.0519],
[0.0319,-3.7043,6370342.0069],
[0.0300,-3.7034,6370342.1510],
[0.0324,-3.7058,6370341.9622],
[0.0297,-3.7052,6370342.0831],
[0.0269,-3.7017,6370342.0821],
[0.0296,-3.7068,6370342.0602],
[0.0282,-3.7045,6370342.0595],
[0.0277,-3.7027,6370342.0495],
[0.0277,-3.7065,6370341.9441],
[0.0260,-3.7048,6370341.9156],
[0.0806,-3.5879,6370341.0791],
[0.0494,-3.6895,6370341.9108],
[-0.0432,-3.8292,6370343.4011],
[-0.1848,-3.6983,6370338.3733],
[0.0594,-3.8016,6370342.2932],
[0.0608,-3.6702,6370341.9858],
[0.0355,-3.7071,6370341.7749],
[0.0386,-3.7082,6370341.8137],
[0.0399,-3.7081,6370342.1022],
[0.0333,-3.7071,6370342.0361],
[0.0391,-3.7095,6370341.9173],
[0.0358,-3.7094,6370341.8443],
[0.0365,-3.7113,6370341.9972],
[0.0342,-3.7109,6370341.9011],
[0.0346,-3.7124,6370341.8879],
[0.0409,-3.7120,6370342.6500],
[0.0336,-3.7135,6370341.8663],
[0.0380,-3.7132,6370342.2314],
[0.0350,-3.7147,6370341.8073],
[0.0393,-3.7141,6370342.3674],
[0.0402,-3.7174,6370341.8738],
[0.0372,-3.7147,6370342.1976],
[0.0301,-3.7092,6370341.9010],
[0.0289,-3.7085,6370341.9100],
[0.0291,-3.7102,6370342.0186],
[0.0310,-3.7096,6370341.9283],
[0.0304,-3.7116,6370342.0336],
[0.0305,-3.7104,6370341.9223],
[0.0302,-3.7131,6370342.0697],
[0.0279,-3.7092,6370342.0813],
[0.0284,-3.7124,6370342.2268],
[0.0310,-3.7146,6370341.9952],
[0.0725,-3.7047,6370342.2131],
[0.0502,-3.7231,6370342.0334],
[-0.2054,-3.7260,6370337.9221],
[0.0227,-3.5384,6370341.3906],
[0.0046,-3.8803,6370384.2876],
[0.0239,-3.8145,6370342.0963],
[0.0370,-3.7179,6370341.8127],
[0.0339,-3.7164,6370341.8756],
[0.0345,-3.7211,6370342.3377],
[0.0340,-3.7184,6370341.9320],
[0.0298,-3.7247,6370342.1810],
[0.0315,-3.7197,6370341.9233],
[0.0301,-3.7162,6370341.7933],
[0.0291,-3.7146,6370341.9160],
[0.0301,-3.7187,6370341.8363],
[0.0320,-3.7172,6370341.7510],
[0.0325,-3.7218,6370342.0195],
[0.0288,-3.7203,6370342.0042],
[0.0264,-3.7138,6370342.1714],
[0.0307,-3.7248,6370342.1724],
[0.0255,-3.7153,6370341.8766],
[0.0280,-3.7144,6370342.0710],
[0.0253,-3.7180,6370341.9647],
[0.0257,-3.7262,6370342.2569],
[0.0266,-3.7208,6370342.2133],
[0.0267,-3.7195,6370342.0901],
[0.0279,-3.7237,6370342.3039],
[0.0277,-3.7221,6370342.2859],
[0.0271,-3.7256,6370342.1752],
[0.0253,-3.7239,6370342.2764],
[0.0248,-3.7166,6370341.8093],
[0.0236,-3.7145,6370342.1043],
[0.0993,-3.5357,6370341.4138],
[-0.1642,-3.7350,6370338.2284],
[0.0887,-3.6336,6370341.6713],
[-0.0622,-3.5125,6370341.3812],
[0.0338,-3.7364,6370342.2736],
[0.0421,-3.7585,6370342.0665],
[0.0301,-3.6958,6370341.7396],
[0.0301,-3.6910,6370341.8955],
[0.0325,-3.6945,6370341.6467],
[0.0285,-3.6961,6370341.8444],
[0.0278,-3.6915,6370341.7856],
[0.0272,-3.6898,6370341.9346],
[0.0275,-3.6939,6370341.7462],
[0.0264,-3.6927,6370341.7551],
[0.0251,-3.6981,6370341.8095],
[0.0270,-3.6957,6370341.8884],
[0.0239,-3.6894,6370341.8679],
[0.0251,-3.7024,6370341.6943],
[0.0241,-3.6956,6370342.0762],
[0.0239,-3.6938,6370342.1742],
[0.0240,-3.6986,6370341.6487],
[0.0228,-3.6972,6370341.7632],
[0.0226,-3.6995,6370341.7943],
[0.0219,-3.6987,6370341.6598],
[0.0244,-3.7011,6370341.7013],
[0.0242,-3.6998,6370341.7086],
[0.0192,-3.6897,6370342.0144],
[0.0301,-3.6976,6370341.8302],
[0.0209,-3.6921,6370342.3147],
[0.0213,-3.6900,6370341.8281],
[0.0211,-3.6943,6370342.2526],
[0.0197,-3.6931,6370342.4030],
[0.0382,-3.6615,6370342.2648],
[0.0214,-3.7015,6370341.6674],
[0.0292,-3.6798,6370341.9412],
[-0.0895,-3.8277,6370341.5461],
[-0.2002,-3.6241,6370379.8460],
[0.0865,-3.7749,6370342.6211],
[0.0267,-3.7070,6370342.2344],
[0.0257,-3.7064,6370342.0731],
[0.0259,-3.7098,6370342.4324],
[0.0253,-3.7077,6370342.3757],
[0.0239,-3.7050,6370342.0186],
[0.0240,-3.7030,6370341.8440],
[0.0218,-3.7059,6370342.1241],
[0.0225,-3.7053,6370342.0635],
[0.0231,-3.7079,6370342.3687],
[0.0241,-3.7059,6370342.1547],
[0.0229,-3.7092,6370342.5419],
[0.0221,-3.7087,6370342.4435],
[0.0228,-3.7103,6370342.5107],
[0.0219,-3.7097,6370342.4810],
[0.0231,-3.7117,6370342.2857],
[0.0220,-3.7112,6370342.3358],
[0.0219,-3.7125,6370342.0797],
[0.0265,-3.7113,6370342.3875],
[0.0214,-3.7037,6370341.7466],
[0.0209,-3.7029,6370341.7326],
[0.0208,-3.7052,6370341.9825],
[0.0195,-3.7040,6370341.8209],
[0.0193,-3.7074,6370342.0706],
[0.0208,-3.7123,6370342.1351],
[0.0202,-3.7108,6370342.1892],
[0.0201,-3.7094,6370342.1157],
[0.0186,-3.7124,6370342.1764],
[0.0176,-3.7068,6370341.9993],
[0.0878,-3.6070,6370341.5292],
[0.0184,-3.7109,6370342.1819],
[0.0124,-3.7853,6370342.3609],
[-0.1640,-3.6048,6370337.4113],
[0.0228,-3.7185,6370342.3522],
[0.0227,-3.7235,6370342.2591],
[0.0245,-3.7219,6370342.2920],
[0.0230,-3.7216,6370342.3073],
[0.0183,-3.7215,6370342.2354],
[0.0222,-3.7262,6370342.2858],
[0.0203,-3.7149,6370342.1801],
[0.0179,-3.7130,6370342.3082],
[0.0213,-3.7169,6370342.2951],
[0.0205,-3.7165,6370342.3274],
[0.0205,-3.7187,6370342.3869],
[0.0201,-3.7179,6370342.4518],
[0.0199,-3.7243,6370342.1636],
[0.0213,-3.7232,6370342.1539],
[0.0171,-3.7159,6370342.3671],
[0.0192,-3.7267,6370342.2858],
[0.0183,-3.7173,6370342.4229],
[0.0185,-3.7161,6370342.3120],
[0.0167,-3.7196,6370342.3213],
[0.0155,-3.7177,6370342.3389],
[0.0182,-3.7229,6370342.2197],
[0.0164,-3.7229,6370342.2684],
[0.0168,-3.7242,6370342.2922],
[0.0170,-3.7253,6370342.4012],
[0.0134,-3.7242,6370342.4318],
[0.0140,-3.7213,6370342.1322],
[0.0768,-3.5719,6370341.4008],
[-0.1070,-3.5225,6370341.0946],
[0.0150,-3.7387,6370342.5029],
[0.0951,-3.6794,6370342.0247],
[0.0172,-3.7638,6370342.1874],
[-0.1911,-3.7746,6370338.1454],
[0.0189,-3.6973,6370341.9159],
[0.0191,-3.6940,6370342.3980],
[0.0185,-3.6983,6370341.9321],
[0.0195,-3.6982,6370341.8108],
[0.0200,-3.7018,6370341.7631],
[0.0202,-3.7006,6370341.7221],
[0.0155,-3.6895,6370341.8599],
[0.0187,-3.7024,6370341.8533],
[0.0152,-3.6942,6370342.6021],
[0.0179,-3.6928,6370342.5879],
[0.0156,-3.6983,6370342.2734],
[0.0183,-3.6960,6370342.2226],
[0.0166,-3.7017,6370342.3008],
[0.0183,-3.7001,6370341.9533],
[0.0121,-3.6915,6370341.9645],
[0.0072,-3.6972,6370342.8944],
[0.0109,-3.6930,6370342.1912],
[0.0138,-3.6924,6370342.0898],
[0.0118,-3.6968,6370342.6058],
[0.0139,-3.6963,6370342.4774],
[0.0127,-3.6987,6370342.4899],
[0.0106,-3.6983,6370342.6425],
[0.0143,-3.7008,6370342.4160],
[0.0144,-3.6993,6370342.3502],
[0.0101,-3.6975,6370342.7290],
[0.0083,-3.6941,6370342.5695],
[-0.1308,-3.5988,6370341.6382],
[0.0078,-3.6609,6370342.0635],
[0.0114,-3.6802,6370342.0595],
[-0.1074,-3.7902,6370339.6557],
[0.0636,-3.7562,6370341.8671],
[0.0688,-3.5458,6370341.1586],
[0.0151,-3.7045,6370342.4307],
[0.0067,-3.7067,6370342.6790],
[0.0152,-3.7071,6370342.4445],
[0.0174,-3.7052,6370342.0015],
[0.0166,-3.7032,6370342.2774],
[0.0163,-3.7076,6370342.2434],
[0.0131,-3.7030,6370342.5742],
[0.0103,-3.7018,6370342.6338],
[0.0131,-3.7052,6370342.5227],
[0.0110,-3.7044,6370342.6691],
[0.0129,-3.7067,6370342.6683],
[0.0147,-3.7059,6370342.4253],
[0.0135,-3.7086,6370342.6116],
[0.0128,-3.7077,6370342.6548],
[0.0075,-3.6999,6370342.8398],
[0.0062,-3.7014,6370342.7758],
[0.0053,-3.7028,6370342.7004],
[0.0082,-3.7009,6370342.7384],
[0.0101,-3.7052,6370342.5921],
[0.0059,-3.7045,6370342.7697],
[0.0077,-3.7079,6370342.5188],
[0.0094,-3.7071,6370342.4356],
[0.0103,-3.7090,6370342.4180],
[0.0048,-3.7086,6370342.8714],
[0.0043,-3.7111,6370342.8738],
[0.0099,-3.7100,6370342.4255],
[0.0349,-3.5051,6370341.8033],
[-0.0173,-3.6783,6370342.4555],
[-0.0039,-3.6914,6370342.7126],
[-0.1369,-3.7699,6370339.4932],
[0.0391,-3.7850,6370341.8158],
[-0.1047,-3.5667,6370341.9005],
[0.0162,-3.7099,6370342.3643],
[0.0069,-3.7161,6370342.5319],
[0.0155,-3.7126,6370342.3307],
[0.0166,-3.7112,6370342.3667],
[0.0161,-3.7159,6370342.3580],
[0.0169,-3.7091,6370342.1898],
[0.0143,-3.7116,6370342.4636],
[0.0132,-3.7105,6370342.4342],
[0.0124,-3.7127,6370342.3731],
[0.0111,-3.7118,6370342.3473],
[0.0062,-3.7123,6370342.7858],
[0.0146,-3.7129,6370342.4554],
[0.0148,-3.7160,6370342.3600],
[0.0128,-3.7155,6370342.2402],
[0.0120,-3.7211,6370342.2277],
[0.0114,-3.7162,6370342.3038],
[0.0086,-3.7114,6370342.4994],
[0.0106,-3.7227,6370342.3558],
[0.0093,-3.7135,6370342.4260],
[0.0070,-3.7118,6370342.6465],
[0.0057,-3.7142,6370342.7178],
[0.0037,-3.7139,6370342.8637],
[0.0068,-3.7180,6370342.5644],
[0.0080,-3.7162,6370342.4026],
[0.0087,-3.7206,6370342.3691],
[0.0092,-3.7183,6370342.2000],
[0.0455,-3.5301,6370341.4567],
[-0.0185,-3.7364,6370342.6227],
[0.0956,-3.7480,6370342.2148],
[-0.1250,-3.5498,6370341.4788],
[-0.1569,-3.7557,6370338.3569],
[-0.0042,-3.7214,6370342.7631]])


antenna_names=['LWA266',
                    'LWA259',
                    'LWA268',
                    'LWA267',
                    'LWA271',
                    'LWA269',
                    'LWA276',
                    'LWA273',
                    'LWA278',
                    'LWA277',
                    'LWA282',
                    'LWA281',
                    'LWA307',
                    'LWA285',
                    'LWA309',
                    'LWA308',
                    'LWA311',
                    'LWA310',
                    'LWA313',
                    'LWA312',
                    'LWA321',
                    'LWA314',
                    'LWA330',
                    'LWA327',
                    'LWA338',
                    'LWA332',
                    'LWA340',
                    'LWA339',
                    'LWA352',
                    'LWA341',
                    'LWA362',
                    'LWA353',
                    'LWA257',
                    'LWA255',
                    'LWA260',
                    'LWA258',
                    'LWA265',
                    'LWA263',
                    'LWA272',
                    'LWA270',
                    'LWA283',
                    'LWA280',
                    'LWA288',
                    'LWA284',
                    'LWA292',
                    'LWA291',
                    'LWA296',
                    'LWA295',
                    'LWA301',
                    'LWA298',
                    'LWA305',
                    'LWA303',
                    'LWA317',
                    'LWA306',
                    'LWA320',
                    'LWA318',
                    'LWA336',
                    'LWA335',
                    'LWA343',
                    'LWA337',
                    'LWA351',
                    'LWA344',
                    'LWA360',
                    'LWA354',
                    'LWA002',
                    'LWA001',
                    'LWA004',
                    'LWA003',
                    'LWA006',
                    'LWA005',
                    'LWA009',
                    'LWA007',
                    'LWA011',
                    'LWA010',
                    'LWA012',
                    'LWA008',
                    'LWA040',
                    'LWA038',
                    'LWA042',
                    'LWA041',
                    'LWA044',
                    'LWA043',
                    'LWA046',
                    'LWA045',
                    'LWA071',
                    'LWA047',
                    'LWA074',
                    'LWA073',
                    'LWA077',
                    'LWA075',
                    'LWA275',
                    'LWA274',
                    'LWA302',
                    'LWA286',
                    'LWA363',
                    'LWA323',
                    'LWA013',
                    'LWA016',
                    'LWA015',
                    'LWA014',
                    'LWA018',
                    'LWA017',
                    'LWA020',
                    'LWA019',
                    'LWA022',
                    'LWA021',
                    'LWA024',
                    'LWA023',
                    'LWA026',
                    'LWA025',
                    'LWA029',
                    'LWA027',
                    'LWA049',
                    'LWA048',
                    'LWA051',
                    'LWA050',
                    'LWA053',
                    'LWA052',
                    'LWA054',
                    'LWA080',
                    'LWA084',
                    'LWA055',
                    'LWA324',
                    'LWA262',
                    'LWA348',
                    'LWA331',
                    'LWA365',
                    'LWA364',
                    'LWA030',
                    'LWA028',
                    'LWA032',
                    'LWA031',
                    'LWA063',
                    'LWA060',
                    'LWA057',
                    'LWA056',
                    'LWA059',
                    'LWA058',
                    'LWA062',
                    'LWA061',
                    'LWA085',
                    'LWA064',
                    'LWA087',
                    'LWA086',
                    'LWA089',
                    'LWA096',
                    'LWA091',
                    'LWA090',
                    'LWA093',
                    'LWA092',
                    'LWA095',
                    'LWA094',
                    'LWA122',
                    'LWA121',
                    'LWA325',
                    'LWA316',
                    'LWA334',
                    'LWA328',
                    'LWA361',
                    'LWA358',
                    'LWA035',
                    'LWA033',
                    'LWA034',
                    'LWA036',
                    'LWA066',
                    'LWA065',
                    'LWA068',
                    'LWA067',
                    'LWA070',
                    'LWA069',
                    'LWA097',
                    'LWA072',
                    'LWA099',
                    'LWA098',
                    'LWA101',
                    'LWA100',
                    'LWA103',
                    'LWA102',
                    'LWA105',
                    'LWA104',
                    'LWA129',
                    'LWA037',
                    'LWA131',
                    'LWA130',
                    'LWA134',
                    'LWA132',
                    'LWA252',
                    'LWA139',
                    'LWA294',
                    'LWA289',
                    'LWA319',
                    'LWA299',
                    'LWA078',
                    'LWA076',
                    'LWA081',
                    'LWA079',
                    'LWA108',
                    'LWA107',
                    'LWA110',
                    'LWA109',
                    'LWA112',
                    'LWA111',
                    'LWA114',
                    'LWA113',
                    'LWA116',
                    'LWA115',
                    'LWA118',
                    'LWA117',
                    'LWA120',
                    'LWA082',
                    'LWA143',
                    'LWA142',
                    'LWA145',
                    'LWA144',
                    'LWA147',
                    'LWA150',
                    'LWA149',
                    'LWA148',
                    'LWA151',
                    'LWA172',
                    'LWA279',
                    'LWA178',
                    'LWA355',
                    'LWA287',
                    'LWA124',
                    'LWA127',
                    'LWA126',
                    'LWA125',
                    'LWA188',
                    'LWA128',
                    'LWA153',
                    'LWA181',
                    'LWA155',
                    'LWA154',
                    'LWA157',
                    'LWA156',
                    'LWA159',
                    'LWA158',
                    'LWA182',
                    'LWA160',
                    'LWA185',
                    'LWA184',
                    'LWA187',
                    'LWA186',
                    'LWA190',
                    'LWA189',
                    'LWA191',
                    'LWA192',
                    'LWA224',
                    'LWA222',
                    'LWA326',
                    'LWA322',
                    'LWA345',
                    'LWA333',
                    'LWA359',
                    'LWA347',
                    'LWA135',
                    'LWA133',
                    'LWA137',
                    'LWA136',
                    'LWA140',
                    'LWA138',
                    'LWA161',
                    'LWA141',
                    'LWA163',
                    'LWA162',
                    'LWA165',
                    'LWA164',
                    'LWA167',
                    'LWA166',
                    'LWA193',
                    'LWA198',
                    'LWA195',
                    'LWA194',
                    'LWA197',
                    'LWA196',
                    'LWA200',
                    'LWA199',
                    'LWA202',
                    'LWA201',
                    'LWA227',
                    'LWA225',
                    'LWA349',
                    'LWA253',
                    'LWA293',
                    'LWA290',
                    'LWA357',
                    'LWA329',
                    'LWA170',
                    'LWA234',
                    'LWA173',
                    'LWA171',
                    'LWA169',
                    'LWA175',
                    'LWA204',
                    'LWA203',
                    'LWA206',
                    'LWA205',
                    'LWA208',
                    'LWA207',
                    'LWA210',
                    'LWA209',
                    'LWA228',
                    'LWA230',
                    'LWA226',
                    'LWA229',
                    'LWA233',
                    'LWA232',
                    'LWA236',
                    'LWA235',
                    'LWA238',
                    'LWA237',
                    'LWA240',
                    'LWA239',
                    'LWA297',
                    'LWA254',
                    'LWA346',
                    'LWA342',
                    'LWA356',
                    'LWA350',
                    'LWA177',
                    'LWA247',
                    'LWA180',
                    'LWA179',
                    'LWA183',
                    'LWA176',
                    'LWA212',
                    'LWA211',
                    'LWA214',
                    'LWA213',
                    'LWA243',
                    'LWA215',
                    'LWA218',
                    'LWA217',
                    'LWA221',
                    'LWA219',
                    'LWA241',
                    'LWA223',
                    'LWA244',
                    'LWA242',
                    'LWA246',
                    'LWA245',
                    'LWA249',
                    'LWA248',
                    'LWA251',
                    'LWA250',
                    'LWA261',
                    'LWA256',
                    'LWA300',
                    'LWA264',
                    'LWA315',
                    'LWA304']