import argparse
from multiprocessing.sharedctypes import Value
import os
from pickle import FALSE
from numpy import isin
import torch
import time
from exp.exp_informer import Exp_Informer
from exp.fc_classify import FC_Classify
from exp.exp_informer_classify import Exp_Informer_Classify
from exp.exp_informer_classify_no_predict import Exp_Informer_Classify_No_Predict
from utils.random_seed import seed_torch
import random
import numpy as np
from utils.log_utils import get_logger
LOGGER = get_logger("informer")


roi_begin = 0
roi_end = 1
cross_validation = 10
itr = 5
for ii in range(0, itr):
    for roi in range(roi_begin, roi_end):
        for cross in range(0,cross_validation):
            cross_s = cross
            model_begin = time.time()
            parser = argparse.ArgumentParser(
                description='[Informer] Long Sequences Forecasting')
            parser.add_argument('--model', type=str, default='informer',
                                help='model of experiment, options: [informer, informerstack, informerlight(TBD)]')
            parser.add_argument('--model_fine_tune', type=str, default='',
                                help='model of experiment, options: [informer, informerstack, informerlight(TBD)]')
            parser.add_argument('--data', type=str, default='SWUWCF', help='data')
            parser.add_argument('--data_test', type=str,
                                default='', help='data for test')
            parser.add_argument('--data_fine_tune', type=str,
                                default='', help='data for fine_tune')
            parser.add_argument('--retest_data_fine_tune', type=str,
                                default='', help='data for fine_tune')

            parser.add_argument('--root_path', type=str,
                                default='./data/ETT/', help='root path of the data file')
            parser.add_argument('--data_path', type=str,
                                default='ETTh1.csv', help='data file')
            parser.add_argument('--features', type=str, default='M',
                                help='forecasting task, options:[M, S, MS]; M:multivariate predict multivariate, S:univariate predict univariate, MS:multivariate predict univariate')  
            parser.add_argument('--target', type=str, default='OT',
                                help='target feature in S or MS task')
            parser.add_argument('--freq', type=str, default='t',
                                help='freq for time features encoding, options:'
                                '[s:secondly, t:minutely, h:hourly, d:daily, b:business days, w:weekly, m:monthly], '
                                'you can also use more detailed freq like 15min or 3h')
            parser.add_argument('--checkpoints', type=str,
                                default='./checkpoints/', help='location of model checkpoints')

            parser.add_argument('--seq_len', type=int, default=64,
                                help='input sequence length of Informer encoder')
            parser.add_argument('--label_len', type=int, default=32,
                                help='start token length of Informer decoder')
            parser.add_argument('--pred_len', type=int, default=8,
                                help='prediction sequence length')  
            parser.add_argument('--enc_in', type=int, default=7,
                                help='encoder input size')
            parser.add_argument('--dec_in', type=int, default=7,
                                help='decoder input size')
            parser.add_argument('--c_out', type=int, default=7, help='output size')

            parser.add_argument('--d_model', type=int,
                                default=454, help='dimension of model')
            parser.add_argument('--l_model', type=int,
                                default=16, help='dimension of model spatial')
            parser.add_argument('--n_heads', type=int,
                                default=8, help='num of heads')
            parser.add_argument('--e_layers', type=int,
                                default=2, help='num of encoder layers')
            parser.add_argument('--d_layers', type=int,
                                default=1, help='num of decoder layers')
            parser.add_argument('--s_layers', type=str,
                                default='3,2,1', help='num of stack encoder layers')
            parser.add_argument('--d_ff', type=int, default=2048,
                                help='dimension of fcn')
            parser.add_argument('--l_ff', type=int, default=1024,
                                help='dimension of fcn spatial')
            parser.add_argument('--factor', type=int, default=5,
                                help='probsparse attn factor')
            parser.add_argument('--factor_spatial', type=int, default=10,
                                help='probsparse attn factor')
            parser.add_argument('--padding', type=int,
                                default=0, help='padding type')
            parser.add_argument('--distil', action='store_false',
                                help='whether to use distilling in encoder, using this argument means not using distilling', default=True)
            parser.add_argument('--dropout', type=float,
                                default=0.05, help='dropout')
            parser.add_argument('--attn', type=str, default='prob',
                                help='attention used in encoder, options:[prob, full]')
            parser.add_argument('--embed', type=str, default='timeF',
                                help='time features encoding, options:[timeF, fixed, learned]')
            parser.add_argument('--activation', type=str,
                                default='gelu', help='activation')
            parser.add_argument('--output_attention', action='store_true',
                                help='whether to output attention in encoder')
            parser.add_argument('--mix', action='store_false',
                                help='use mix attention in generative decoder', default=True)
            parser.add_argument('--cols', type=str, nargs='+',
                                help='certain cols from the data files as the input features')
            parser.add_argument('--num_workers', type=int,
                                default=0, help='data loader num workers')
            parser.add_argument('--itr', type=int, default=1,
                                help='experiments times')
            parser.add_argument('--train_epochs', type=int,
                                default=20, help='train epochs')
            parser.add_argument('--batch_size', type=int, default=128,
                                help='batch size of train input data')  # 改 默认32
            parser.add_argument('--patience', type=int, default=3,
                                help='early stopping patience')
            parser.add_argument('--learning_rate', type=float,
                                default=0.0001, help='optimizer learning rate')
            parser.add_argument('--des', type=str,
                                default='test', help='exp description')
            parser.add_argument('--loss', type=str,
                                default='mse', help='loss function')
            parser.add_argument('--lradj', type=str,
                                default='type1', help='adjust learning rate')
            parser.add_argument('--use_amp', action='store_true',
                                help='use automatic mixed precision training', default=False)
            parser.add_argument('--inverse', action='store_true',
                                help='inverse output data', default=False)  
            parser.add_argument('--test_part', type=str,
                                default='test', help='test part, test train valid all')                    
            parser.add_argument('--do_train', action='store_true',
                                help='whether to train model', default=False)
            parser.add_argument('--do_test', action='store_true',
                                help='whether to train model', default=False)
            parser.add_argument('--do_retest', action='store_true',
                                help='whether to train model', default=False)
            parser.add_argument('--do_predict', action='store_true',
                                help='whether to predict unseen future data', default=False)

            parser.add_argument('--use_gpu', action='store_true',
                                default=False, help='use gpu')
            parser.add_argument('--use_multi_gpu', action='store_true',
                                help='use multiple gpus', default=False)
            parser.add_argument('--devices', type=str,
                                default='0,1,2,3', help='device ids of multile gpus')
            parser.add_argument('--class_type', type=str,
                                default='FC_Classify', help='the class of exp')
            parser.add_argument('--ML_model', type=str,
                                default='svm', help='the method of ML model')
            parser.add_argument('--fc_state', type=str,
                                default='dynamic', help='the state of FC dynamic or static')
            parser.add_argument('--class_output', type=str,
                                default=4, help='the output dimension of model')
            parser.add_argument('--checkpoint_name', type=str,
                                default="checkpoint.pth", help='the name of checkpoint for other fine-tuning')

            args = parser.parse_args()

            if args.use_gpu:
                args.devices = args.devices.replace(' ', '')
                device_ids = args.devices.split(',')
                if args.use_multi_gpu:
                    args.device_ids = [int(id_) for id_ in device_ids]
                else:
                    if args.devices:
                        args.device_ids = int(device_ids[0])
                    else:
                        raise ValueError('Please input device id!')

            data_parser = {
                'SWUWCF': {'data': 'SWUWCF.hdf5', 'T': 0, 'M': [454, 454, 454], 'S': [1, 1, 1], 'MS': [454, 454, 1]},
                'MDD_4task': {'data': 'MDD_4task.hdf5', 'T': 0, 'M': [454, 454, 454], 'S': [1, 1, 1], 'MS': [454, 454, 1]},
                'MULTI_RUMI_EMO': {'data': 'MULTI_RUMI_EMO.hdf5', 'T': 0, 'M': [454, 454, 454], 'S': [1, 1, 1], 'MS': [454, 454, 1]},
                'Rumination_Suzhou': {'data': 'Rumination_Suzhou.hdf5', 'T': 0, 'M': [454, 454, 454], 'S': [1, 1, 1], 'MS': [454, 454, 1]},
            }       
                        
            if args.data_fine_tune and args.data_test:
                LOGGER.error("Error: data_fine_tune and data_test can not be used at the same time!")
            else:
                if args.data_fine_tune:
                    use_data = args.data_fine_tune
                elif args.data_test:
                    use_data = args.data_test
                else:
                    use_data = args.data
            if use_data in data_parser.keys():
                data_info = data_parser[use_data]
                args.data_path = data_info['data']
                LOGGER.info("data path:{}".format(args.data_path))
                args.target = roi
                args.enc_in, args.dec_in, args.c_out = data_info[args.features]
            else:
                raise ValueError("Error: data not found!")

            args.s_layers = [int(s_l)
                            for s_l in args.s_layers.replace(' ', '').split(',')]
            args.detail_freq = args.freq
            args.freq = args.freq[-1:]
                
            LOGGER.info('Args in experiment:')
            LOGGER.info(args)
            Class_dict = {
                'Exp_Informer':Exp_Informer,
                'Exp_Informer_Classify':Exp_Informer_Classify,
            }
            Exp = Class_dict[args.class_type]

            LOGGER.info("cross: {}  itr: {}".format(cross, ii))
            args.setting = '{}_{}_ft{}_sl{}_ll{}_pl{}_dm{}_nh{}_el{}_dl{}_df{}_at{}_fc{}_eb{}_dt{}_mx{}_{}_{}_{}'.format(args.model, args.data, args.features,
                                                                                                                        args.seq_len, args.label_len, args.pred_len,
                                                                                                        args.d_model, args.n_heads, args.e_layers, args.d_layers, args.d_ff, args.attn, args.factor,
                                                                                                                        args.embed, args.distil, args.mix, args.des, ii) 
            exp = Exp(args)
            if args.do_train:
                LOGGER.info('>>>>>>>start training : {}>>>>>>>>>>>>>>>>>>>>>>>>>>'.format(
                    args.setting))
                exp.train(cross,ii)

            if args.do_test:
                LOGGER.info('>>>>>>>testing : {}<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'.format(
                    args.setting))
                exp.test(cross,ii)

            if args.do_retest:
                LOGGER.info('>>>>>>>retesting : {}<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'.format(
                    args.setting))
                exp.retest(cross, ii, flag=args.test_part)

            if args.do_predict:
                LOGGER.info('>>>>>>>predicting : {}<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'.format(
                    args.setting))
                exp.predict(cross, ii, True)

            torch.cuda.empty_cache()
            model_end = time.time()
            LOGGER.info("time:{}".format(model_end-model_begin))
