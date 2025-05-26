from sklearn.metrics import confusion_matrix
from data.data_loader import Dataset_IPCASS3, Dataset_Multi_Timepoint, Dataset_Multi_Task
from exp.exp_basic import Exp_Basic
from models.model import Informer, InformerStack , InformerClassify, InformerSpatial, InformerSpatialClassify

from utils.tools import EarlyStopping, adjust_learning_rate
from utils.metrics import MAE, MSE, RMSE, RSE, metric, CORR, ACC, PRE, RECALL, F1, CONFUSION_MATRIX, CLASS_REPORT

import numpy as np
import pandas as pd

import torch
import torch.nn as nn
from torch import optim
from torch.utils.data import DataLoader
from torch.utils.tensorboard import SummaryWriter 

from tqdm import tqdm
import torch.nn.functional as F
from utils.log_utils import get_logger
import h5py




import os
import time

import warnings
LOGGER = get_logger("informer")

class Exp_Informer(Exp_Basic):
    def __init__(self, args):
        super(Exp_Informer, self).__init__(args)
    def _load_model(self, best_model_path):
        if self.args.use_gpu:
            state_dict = torch.load(best_model_path)
        else:
            state_dict = torch.load(best_model_path, map_location=torch.device('cpu'))
        if not self.args.use_multi_gpu:
            from collections import OrderedDict
            new_state_dict = OrderedDict()
            for k, v in state_dict.items():
                name = k.replace('.module.','.') 
                name = k.replace('module.','') 
                new_state_dict[name] = v
            state_dict = new_state_dict
        return state_dict
    def _build_model(self):
        model_dict = {
            'informer':Informer,
            'informer_stack':InformerStack,
            'informer_spatial':InformerSpatial,
        }
        if self.args.model_fine_tune:
            model_use = self.args.model_fine_tune
        else:
            model_use = self.args.model
        if model_use=='informer' or model_use=='informer_stack':
            e_layers = self.args.e_layers if model_use!='informer_stack' else self.args.s_layers
            model = model_dict[model_use](
                self.args.enc_in,
                self.args.dec_in, 
                self.args.c_out, 
                self.args.seq_len, 
                self.args.label_len,
                self.args.pred_len, 
                self.args.factor,
                self.args.d_model, 
                self.args.n_heads, 
                e_layers,
                self.args.d_layers, 
                self.args.d_ff,
                self.args.l_model,
                self.args.l_ff,
                self.args.dropout, 
                self.args.attn,
                self.args.embed,
                self.args.freq,
                self.args.activation,
                self.args.output_attention,
                self.args.distil,
                self.args.mix,
                self.device
            ).float()
        elif model_use=='informer_spatial':
            model = model_dict[model_use](
                self.args.enc_in,
                self.args.dec_in, 
                self.args.c_out, 
                self.args.seq_len, 
                self.args.label_len,
                self.args.pred_len, 
                self.args.factor,
                self.args.factor_spatial,
                self.args.d_model, 
                self.args.n_heads, 
                self.args.e_layers, 
                self.args.d_layers, 
                self.args.d_ff,
                self.args.l_model,
                self.args.l_ff,
                self.args.dropout, 
                self.args.attn,
                self.args.embed,
                self.args.freq,
                self.args.activation,
                self.args.output_attention,
                self.args.distil,
                self.args.mix,
                self.device
            ).float()

        if self.args.use_multi_gpu and self.args.use_gpu:
            LOGGER.info("{}".format(self.args.device_ids))
            model = nn.DataParallel(model)
        if self.args.data_fine_tune: 
            path = os.path.join(self.args.checkpoints, self.args.setting)
            best_model_path = path+'/'+'checkpoint.pth'
            save_model = self._load_model(best_model_path)
            model_dict =  model.state_dict()
            state_dict = {k:v for k,v in save_model.items() if k in model_dict.keys()}
            model_dict.update(state_dict)
            model.load_state_dict(model_dict)
        return model

    def _get_data(self, cross, flag):
        args = self.args

        data_dict = {
            'SWUWCF':Dataset_Multi_Timepoint,
            "MULTI_RUMI_EMO":Dataset_Multi_Timepoint,
            "Rumination_Suzhou":Dataset_Multi_Timepoint,
        }

        if self.args.data_fine_tune and self.args.data_test:
            raise ValueError("Error: data_fine_tune and data_test can not be used at the same time!")

        else:
            if self.args.data_fine_tune:
                Data = data_dict[self.args.data_fine_tune]
            elif self.args.data_test:
                Data = data_dict[self.args.data_test]
            else:
                Data = data_dict[self.args.data]
        timeenc = 0 if args.embed!='timeF' else 1

        if flag == 'train':
            shuffle_flag = True; drop_last = True; batch_size = args.batch_size; freq=args.freq
        elif flag == 'val':
            shuffle_flag = True; drop_last = True; batch_size = args.batch_size; freq=args.freq
        elif flag=='pred':
            shuffle_flag = False; drop_last = False; batch_size = 1; freq=args.detail_freq
            Data = Dataset_Pred
        else: 
            shuffle_flag = False; drop_last = True; batch_size = args.batch_size; freq=args.freq 

        data_set = Data(
            root_path=args.root_path,
            data_path=args.data_path,
            flag=flag,
            size=[args.seq_len, args.label_len, args.pred_len],
            features=args.features,
            target=args.target,
            inverse=args.inverse,
            timeenc=timeenc,
            freq=freq,
            cols=args.cols,
            data_fine_tune = args.data_fine_tune,
            setting = args.setting,
            cross = cross
        ) 
        LOGGER.info("{},{}".format(flag, len(data_set))) 
        data_loader = DataLoader(
            data_set,
            batch_size=batch_size,
            shuffle=shuffle_flag, 
            num_workers=args.num_workers,
            drop_last=drop_last)

        return data_set, data_loader

    def _select_optimizer(self):
        model_optim = optim.Adam(self.model.parameters(), lr=self.args.learning_rate)
        return model_optim
    
    def _select_criterion_predict(self):
        criterion =  nn.MSELoss()
        return criterion
    def _select_criterion_classify(self):
        criterion =  nn.CrossEntropyLoss()
        return criterion

    def vali(self, vali_data, vali_loader, criterion,flag):
        self.model.eval()
        total_loss = []
        preds = []
        trues = []

        for i, (batch_x,batch_y) in enumerate(vali_loader):
            pred, true = self._process_one_batch(
                vali_data, batch_x, batch_y)
            loss = criterion(pred.detach().cpu(), true.detach().cpu())
            total_loss.append(loss)
            preds.append(pred.detach().cpu().numpy())
            trues.append(true.detach().cpu().numpy())

        preds = np.array(preds)
        trues = np.array(trues)
        preds = preds.reshape(-1, preds.shape[-2], preds.shape[-1])
        trues = trues.reshape(-1, trues.shape[-2], trues.shape[-1])
        corr_sub = np.zeros((vali_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        
        for sub in range(vali_data.data_num):
            corr_board0 = sum(vali_data.num_win[0:sub])
            corr_board1 = sum(vali_data.num_win[0:sub+1])
            if sub != vali_data.data_num-1:
                corr_sub[sub] = CORR(preds[corr_board0:corr_board1],trues[corr_board0:corr_board1])
            else:
                corr_sub[sub] = CORR(preds[corr_board0:],trues[corr_board0:])

        corr_mean2 = np.mean(corr_sub)

        LOGGER.info("corr_mean:{}".format(corr_mean2))
        corr = CORR(preds, trues)
        total_corr = np.mean(corr)
        total_loss = np.average(total_loss)
        self.model.train()
        return corr_mean2


    def train(self, cross,ii):
        train_data, train_loader = self._get_data(cross, flag = 'train')
        vali_data, vali_loader = self._get_data(cross, flag = 'val')
        torch.set_num_threads(4)

        path = os.path.join(self.args.checkpoints, self.args.setting)
        if not os.path.exists(path):
            os.makedirs(path)
        writer = SummaryWriter('./tensorboard')
        time_now = time.time()
        
        train_steps = len(train_loader)
        early_stopping = EarlyStopping(patience=self.args.patience, verbose=True)
        
        model_optim = self._select_optimizer()
        criterion =  self._select_criterion_predict()

        if self.args.use_amp:
            scaler = torch.cuda.amp.GradScaler()
        tensorboard_num  = 0 
        for epoch in range(self.args.train_epochs):
            iter_count = 0
            train_loss = []
            
            self.model.train()
            epoch_time = time.time()
            for i, (batch_x,batch_y) in enumerate(train_loader): 
                iter_count += 1
                
                model_optim.zero_grad()
                if self.args.output_attention:
                    pred, true,attn = self._process_one_batch(
                        train_data, batch_x, batch_y)
                else:
                    pred, true = self._process_one_batch(
                        train_data, batch_x, batch_y)
                loss = criterion(pred, true)
                train_loss.append(loss.item())
                writer.add_scalar('loss_train', loss.item(), tensorboard_num) 
                tensorboard_num  = tensorboard_num+1 
                if (i+1) % 1000==0:
                    vali_loss = self.vali(vali_data, vali_loader, criterion,flag = "val")
                    writer.add_scalar('loss_val', vali_loss, tensorboard_num)
                    time_now = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                    LOGGER.info("{}".format(time_now))
                    time_nows = time.time()
                    LOGGER.info("\titers: {0}, epoch: {1} | loss: {2:.7f} | valid_corr{3:.7f}".format(i + 1, epoch + 1, loss.item(), vali_loss))
                    speed = (time.time()-time_nows)/iter_count
                    left_time = speed*((self.args.train_epochs - epoch)*train_steps - i)
                    LOGGER.info('\tspeed: {:.4f}s/iter; left time: {:.4f}s'.format(speed, left_time))
                    iter_count = 0
                elif (i+1) % 100==0:
                    time_now = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                    LOGGER.info("{}".format(time_now))
                    time_nows = time.time()
                    LOGGER.info("\titers: {0}, epoch: {1} | loss: {2:.7f}".format(i + 1, epoch + 1, loss.item()))
                    speed = (time.time()-time_nows)/iter_count
                    left_time = speed*((self.args.train_epochs - epoch)*train_steps - i)
                    LOGGER.info('\tspeed: {:.4f}s/iter; left time: {:.4f}s'.format(speed, left_time))
                    iter_count = 0

                if self.args.use_amp:
                    scaler.scale(loss).backward()
                    scaler.step(model_optim)
                    scaler.update()
                else:
                    loss.backward()
                    model_optim.step()

            LOGGER.info("Epoch: {} cost time: {}".format(epoch+1, time.time()-epoch_time))
            train_loss = np.average(train_loss)
            vali_loss = self.vali(vali_data, vali_loader, criterion,flag = "val")
 
            LOGGER.info("Epoch: {0}, Steps: {1} | Train Loss: {2:.7f} Vali Loss: {3:.7f}".format(
                epoch + 1, train_steps, train_loss, vali_loss))
            early_stopping(vali_loss, self.model, path, self.args.data_fine_tune, cross)
            if early_stopping.early_stop:
                LOGGER.info("Early stopping")
                break

            adjust_learning_rate(model_optim, epoch+1, self.args)
        if self.args.data_fine_tune:
            best_model_path = path+'/'+self.args.data_fine_tune+str(cross)+'checkpoint.pth'
        else:
            best_model_path = path+'/'+str(cross)+'checkpoint.pth'
        state_dict = self._load_model(best_model_path)
        self.model.load_state_dict(state_dict)        
        return self.model

    def test(self, cross, ii):
        test_data, test_loader = self._get_data(cross,flag='test')
        
        self.model.eval()
        
        preds = []
        trues = []
        writer = SummaryWriter('./tensorboard')
        tensorboard_num = 0
        for i, (batch_x,batch_y) in enumerate(test_loader):
            pred, true = self._process_one_batch(
                test_data, batch_x, batch_y)
            preds.append(pred.detach().cpu().numpy())
            trues.append(true.detach().cpu().numpy())

        preds = np.array(preds)
        trues = np.array(trues)
        LOGGER.info('test shape:{},{}'.format(preds.shape, trues.shape))
        preds = preds.reshape(-1, preds.shape[-2], preds.shape[-1])
        trues = trues.reshape(-1, trues.shape[-2], trues.shape[-1])
        LOGGER.info('test shape:{},{}'.format(preds.shape, trues.shape))

        folder_path = './results/' + self.args.setting +'/'
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        np.save(folder_path+'pred.npy', preds)
        np.save(folder_path+'true.npy', trues)

        mae, mse, rmse, mape, mspe, rse, corr = metric(preds, trues)
        LOGGER.info('mse:{}, mae:{},rmse{},mape{},mspe{},rse{},corr{}'.format(mse, mae, rmse, mape, mspe, rse, corr))
        np.save(folder_path+'metrics.npy', np.array([mae, mse, rmse, mape, mspe,rse,corr]))
       
        LOGGER.info("本次实验测试结果：{}".format(preds))
        LOGGER.info("真实值：{}".format(trues))

        return

    def predict(self, cross, load=False):
        pred_data, pred_loader = self._get_data(cross,flag='pred')
        
        if load:
            path = os.path.join(self.args.checkpoints, self.args.setting)
            if self.args.data_fine_tune:
                best_model_path = path+'/'+self.args.data_fine_tune+str(cross)+'checkpoint.pth'
            else:
                best_model_path = path+'/'+str(cross)+'checkpoint.pth'
            state_dict = self._load_model(best_model_path)
            self.model.load_state_dict(state_dict)

        self.model.eval()
        
        preds = []
        
        for i, (batch_x,batch_y) in enumerate(tqdm(pred_loader)):
            pred, true = self._process_one_batch(
                pred_data, batch_x, batch_y)
            preds.append(pred.detach().cpu().numpy())

        preds = np.array(preds)
        preds = preds.reshape(-1, preds.shape[-2], preds.shape[-1])
        
        folder_path = './results/' + self.args.setting +'/'
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        np.save(folder_path+'real_prediction.npy', preds)
        LOGGER.info("本次实验预测未来结果：{}".format(preds))

        return

    def _process_one_batch(self, dataset_object, batch_x, batch_y):
        batch_x = batch_x.float().to(self.device)
        batch_y = batch_y.float()

        if self.args.padding==0:
            dec_inp = torch.zeros([batch_y.shape[0], self.args.pred_len, batch_y.shape[-1]]).float()
        elif self.args.padding==1:
            dec_inp = torch.ones([batch_y.shape[0], self.args.pred_len, batch_y.shape[-1]]).float()
        dec_inp = torch.cat([batch_y[:,:self.args.label_len,:], dec_inp], dim=1).float().to(self.device)
        if self.args.use_amp:
            with torch.cuda.amp.autocast():
                if self.args.output_attention:
                    outputs,attn = self.model(batch_x, dec_inp)
                else:
                    outputs = self.model(batch_x, dec_inp)
        else:
            if self.args.output_attention:
                outputs,attn = self.model(batch_x, dec_inp)
            else:
                outputs = self.model(batch_x, dec_inp)
        if self.args.inverse:
            outputs = dataset_object.inverse_transform(outputs)
        f_dim = -1 if self.args.features=='MS' else 0
        batch_y = batch_y[:,-self.args.pred_len:,f_dim:].to(self.device)
        if self.args.output_attention:
            return outputs, batch_y, attn
        else:
            return outputs, batch_y


    def retest(self, cross, ii, flag):

        test_data, test_loader = self._get_data(cross, flag)
        torch.set_num_threads(4)

        path = os.path.join(self.args.checkpoints, self.args.setting)
        if self.args.retest_data_fine_tune:
            best_model_path = path+'/'+self.args.retest_data_fine_tune+str(cross)+'checkpoint.pth'
        else:
            best_model_path = path+'/'+'0'+'checkpoint.pth'
        state_dict = self._load_model(best_model_path)
        self.model.load_state_dict(state_dict)
        self.model.eval()
        
        preds = []
        trues = []
        batch_xs = []
        if self.args.output_attention:
            attns=torch.zeros(sum(test_data.num_win),test_data.seq_x.shape[1],test_data.seq_x.shape[1])
        
        for i, (batch_x,batch_y) in enumerate(tqdm(test_loader)):
            if self.args.output_attention:
                pred, true,attn = self._process_one_batch(
                    test_data, batch_x, batch_y)
            else:
                pred, true = self._process_one_batch(
                    test_data, batch_x, batch_y)
            preds.append(pred.detach().cpu().numpy())
            trues.append(true.detach().cpu().numpy())
            batch_xs.append(batch_x.detach().cpu().numpy())
            if self.args.output_attention:
                if i != int((sum(test_data.num_win)/test_loader.batch_size)//1):
                    attns[i*test_loader.batch_size:(i+1)*test_loader.batch_size]=attn[0].detach().mean(1).cpu()
                else:
                    attns[i*test_loader.batch_size:]=attn[0].detach().mean(1).cpu()
        

        preds = np.array(preds)
        trues = np.array(trues)
        LOGGER.info('test shape:{}, {}'.format(preds.shape, trues.shape))
        preds = preds.reshape(-1, preds.shape[-2], preds.shape[-1])
        trues = trues.reshape(-1, trues.shape[-2], trues.shape[-1])
        LOGGER.info('test shape:{}, {}'.format(preds.shape, trues.shape))
        
        folder_path = './results/attn/' +self.args.setting +'/'
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        corr_sub = np.zeros((test_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        mse_sub = np.zeros((test_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        mae_sub = np.zeros((test_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        rmse_sub = np.zeros((test_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        rse_sub = np.zeros((test_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        if self.args.output_attention:
            attn_sub = torch.zeros((test_data.data_num,attns.shape[1],attns.shape[2]))

        for sub in range(test_data.data_num):
            corr_board0 = sum(test_data.num_win[0:sub])
            corr_board1 = sum(test_data.num_win[0:sub+1])
            if sub != test_data.data_num-1:
                corr_sub[sub] = CORR(preds[corr_board0:corr_board1],trues[corr_board0:corr_board1])
                mse_sub[sub] = MSE(preds[corr_board0:corr_board1],trues[corr_board0:corr_board1])
                mae_sub[sub] = MAE(preds[corr_board0:corr_board1],trues[corr_board0:corr_board1])
                rmse_sub[sub] = RMSE(preds[corr_board0:corr_board1],trues[corr_board0:corr_board1])
                rse_sub[sub] = RSE(preds[corr_board0:corr_board1],trues[corr_board0:corr_board1])
                if self.args.output_attention:
                    attn_sub[sub] = attns[corr_board0:corr_board1].mean(0)
            else:
                corr_sub[sub] = CORR(preds[corr_board0:],trues[corr_board0:])
                mse_sub[sub] = MSE(preds[corr_board0:],trues[corr_board0:])
                mae_sub[sub] = MAE(preds[corr_board0:],trues[corr_board0:])
                rmse_sub[sub] = RMSE(preds[corr_board0:],trues[corr_board0:])
                rse_sub[sub] = RSE(preds[corr_board0:],trues[corr_board0:])
                if self.args.output_attention:
                    attn_sub[sub] = attns[corr_board0:].mean(0)

        corr_mean2 = corr_sub.mean(0).mean(0).mean(0)
        mse_mean2 = mse_sub.mean(0).mean(0).mean(0)
        mae_mean2 = mae_sub.mean(0).mean(0).mean(0)
        rmse_mean2 = mse_sub.mean(0).mean(0).mean(0)
        rse_mean2 = rse_sub.mean(0).mean(0).mean(0)

        LOGGER.info("corr_mean:{}".format(corr_mean2))
        LOGGER.info("mse_mean:{}".format(mse_mean2))
        LOGGER.info("mae_mean:{}".format(mae_mean2))
        LOGGER.info("rmse_mean:{}".format(rmse_mean2))
        LOGGER.info("rse_mean:{}".format(rse_mean2))


        if self.args.retest_data_fine_tune:
            if self.args.output_attention:
                np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'batch_xs.npy', batch_xs)
                np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'attn.npy', attn_sub)
        else:
            np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+str(cross)+'repred.npy', preds)
            np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+str(cross)+'retrue.npy', trues)
            np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+str(cross)+'remetrics_onlypredict.npy', np.array([mse_sub, mae_sub, rmse_sub, rse_sub,corr_sub]))
            np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+str(cross)+'remetrics_onlypredict_mean.npy', np.array([mse_mean2, mae_mean2, rmse_mean2, rse_mean2,corr_mean2]))
            if self.args.output_attention:

                f1 = h5py.File(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+str(cross)+'attn_sub.hdf5', 'w')

                sub_i = 0
                for sub_i in range(attn_sub.shape[0]):
                        f1[str(sub_i)] = attns[sub_i]
                        sub_i = sub_i+1
                f1.close()
        return