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



import os
import time

import warnings
LOGGER = get_logger("informer")

class Exp_Informer_Classify(Exp_Basic):
    def __init__(self, args):
        super(Exp_Informer_Classify, self).__init__(args)
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
            'informer_classify':InformerClassify,
            'informer_spatial_classify':InformerSpatialClassify
        }
        if self.args.model_fine_tune:
            model_use = self.args.model_fine_tune
        else:
            model_use = self.args.model
        if model_use=='informer_classify':
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
                self.device,
                self.args.class_output,
            ).float()
        elif model_use=='informer_spatial_classify':
            e_layers = self.args.e_layers if model_use!='informer_stack' else self.args.s_layers
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
                self.device,
                self.args.class_output,
            ).float()

        if self.args.use_multi_gpu and self.args.use_gpu:
            LOGGER.info("{}".format(self.args.device_ids))
            model = nn.DataParallel(model)
        if self.args.data_fine_tune:
            path = os.path.join(self.args.checkpoints, self.args.setting)
            best_model_path = path+'/'+self.args.checkpoint_name
            save_model = self._load_model(best_model_path)
            model_dict =  model.state_dict()
            if self.args.checkpoint_name == '0checkpoint.pth':
                state_dict = {k:v for k,v in save_model.items() if k in model_dict.keys()}
            else:
                state_dict = {k:v for k,v in save_model.items() if k in model_dict.keys() and k not in ["module.projection_classify.weight","module.projection_classify.bias"]}
            model_dict.update(state_dict)
            model.load_state_dict(model_dict)
        return model

    def _get_data(self, cross, flag):
        args = self.args

        data_dict = {
            "MDD_4task":Dataset_Multi_Task,
            "MULTI_RUMI_EMO":Dataset_Multi_Task,
            "Rumination_Suzhou":Dataset_Multi_Task,
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

    def vali(self, vali_data, vali_loader, criterion_predict,criterion_classify, lamda, flag):
        self.model.eval()
        total_loss = []
        preds = []
        trues = []
        pred_labels = []
        batch_labels = []
        for i, (batch_x, batch_y, batch_label, batch_station) in enumerate(vali_loader):
            if self.args.output_attention:
                    pred, true, pred_label,true_label,attn = self._process_one_batch(
                        vali_data, batch_x, batch_y, batch_label, batch_station)
            else:
                    pred, true, pred_label,true_label = self._process_one_batch(
                        vali_data, batch_x, batch_y, batch_label, batch_station)  

            loss_predict = criterion_predict(pred.detach().cpu(), true.detach().cpu())
            loss_classify = criterion_classify(pred_label.detach().cpu(),batch_label.detach().cpu())
            loss = lamda*loss_predict + (1-lamda)*loss_classify
            
            total_loss.append(loss)
            preds.append(pred.detach().cpu().numpy())
            trues.append(true.detach().cpu().numpy())
            pred_labels.append(pred_label.detach().cpu().numpy())
            batch_labels.append(batch_label.detach().cpu().numpy())

        preds = np.array(preds)
        trues = np.array(trues)
        preds = preds.reshape(-1, preds.shape[-2], preds.shape[-1])
        trues = trues.reshape(-1, trues.shape[-2], trues.shape[-1])
        pred_labels = np.array(pred_labels)
        batch_labels = np.array(batch_labels)
        pred_labels = pred_labels.reshape(-1, pred_labels.shape[-1])
        pred_labels = np.argmax(pred_labels, axis=1)
        pred_labels = pred_labels.reshape(-1, 1)
        batch_labels = batch_labels.reshape(-1, 1)
        corr_sub = np.zeros((vali_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        acc_sub = np.zeros((vali_data.data_num,pred_labels.shape[1]),dtype = None, order = 'C')

        
        for sub in range(vali_data.data_num):
            corr_board0 = sum(vali_data.num_win[0:sub])
            corr_board1 = sum(vali_data.num_win[0:sub+1])
            if sub != vali_data.data_num-1:
                corr_sub[sub] = CORR(preds[corr_board0:corr_board1],trues[corr_board0:corr_board1])
                acc_sub[sub] = ACC(pred_labels[corr_board0:corr_board1],batch_labels[corr_board0:corr_board1])

            else:
                corr_sub[sub] = CORR(preds[corr_board0:],trues[corr_board0:])
                acc_sub[sub] = ACC(pred_labels[corr_board0:],batch_labels[corr_board0:])

        corr_mean2 = np.mean(corr_sub)
        acc_mean2 = np.mean(acc_sub)
        LOGGER.info("corr_mean:{} | acc_mean:{}".format(corr_mean2,acc_mean2))
        self.model.train()
        return corr_mean2, acc_mean2


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
        criterion_predict =  self._select_criterion_predict()
        criterion_classify =  self._select_criterion_classify()


        if self.args.use_amp:
            scaler = torch.cuda.amp.GradScaler()
        tensorboard_num  = 0 
        for epoch in range(self.args.train_epochs):
            iter_count = 0
            train_loss = []
            train_loss_predict = []
            train_loss_classify = []
            
            self.model.train()
            epoch_time = time.time()
            time_nows = time.time()

            for i, (batch_x, batch_y, batch_label, batch_station) in enumerate(train_loader): 

                model_optim.zero_grad()
                if self.args.output_attention:
                    pred, true, pred_label,true_label,attn = self._process_one_batch(
                        train_data, batch_x, batch_y, batch_label, batch_station)
                else:
                    pred, true, pred_label,true_label = self._process_one_batch(
                        train_data, batch_x, batch_y, batch_label, batch_station)  
                loss_predict = criterion_predict(pred, true)
                loss_classify = criterion_classify(pred_label,batch_label.to(self.device))
                lamda = 0.6
                loss = lamda*loss_predict + (1-lamda)*loss_classify



                train_loss.append(loss.item())
                train_loss_predict.append(loss_predict.item())
                train_loss_classify.append(loss_classify.item())

                writer.add_scalar('loss_train', loss.item(), tensorboard_num) 
                writer.add_scalar('loss_train_predict', loss_predict.item(), tensorboard_num) 
                writer.add_scalar('loss_train_classify', loss_classify.item(), tensorboard_num) 

                tensorboard_num = tensorboard_num+1
                if (i+1) % 5000==0:
                    vali_corr, vali_acc = self.vali(vali_data, vali_loader, criterion_predict, criterion_classify,lamda, flag = "val")
                    writer.add_scalar('vali_corr', vali_corr, tensorboard_num) 
                    writer.add_scalar('vali_acc', vali_acc, tensorboard_num) 
                    time_now = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                    LOGGER.info("{}".format(time_now))
                    LOGGER.info("\titers: {0}, epoch: {1} | loss: {2:.7f} val_corr: {3:.7f} val_acc: {4:.7f}".format(i + 1, epoch + 1, loss.item(), vali_corr.item(), vali_acc.item()))
                    LOGGER.info("\tloss_prdict: {0:.7f} loss_classify: {1:.7f}".format(loss_predict.item(), loss_classify.item()))
                    
                    speed = (time.time()-time_nows)/i
                    left_time = speed*((self.args.train_epochs - epoch)*train_steps - i)
                    LOGGER.info('\tspeed: {:.4f}s/iter; left time: {:.4f}s'.format(speed, left_time))
                elif (i+1) % 100==0:
                    time_now = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                    LOGGER.info("{}".format(time_now))
                    time_nows = time.time()
                    LOGGER.info("\titers: {0}, epoch: {1} | loss: {2:.7f}".format(i + 1, epoch + 1, loss.item()))
                    LOGGER.info("\tloss_prdict: {0:.7f} loss_classify: {1:.7f}".format(loss_predict.item(), loss_classify.item()))

                    speed = (time.time()-time_nows)/i
                    left_time = speed*((self.args.train_epochs - epoch)*train_steps - i)
                    LOGGER.info('\tspeed: {:.4f}s/iter; left time: {:.4f}s'.format(speed, left_time))
                iter_count = iter_count


                
                if self.args.use_amp:
                    scaler.scale(loss).backward()
                    scaler.step(model_optim)
                    scaler.update()
                else:
                    loss.backward()
                    model_optim.step()

            LOGGER.info("Epoch: {} cost time: {}".format(epoch+1, time.time()-epoch_time))
            train_loss = np.average(train_loss)
            train_loss_predict = np.average(train_loss_predict)
            train_loss_classify = np.average(train_loss_classify)
            vali_corr, vali_acc = self.vali(vali_data, vali_loader, criterion_predict, criterion_classify,lamda, flag = "val")
            LOGGER.info("Epoch: {0}, Steps: {1} | Train Loss: {2:.7f} Train Lodss predict: {3:.7f} Train Loss classify: {4:.7f} Vali Corr: {5:.7f} Vali Acc: {6:.7f}".format(
                epoch + 1, train_steps, train_loss, train_loss_predict, train_loss_classify, vali_corr, vali_acc))
            early_stop_index = lamda*vali_corr + (1-lamda)*vali_acc

            early_stopping(early_stop_index, self.model, path, self.args.data_fine_tune, cross)
            if early_stopping.early_stop:
                LOGGER.info("Early stopping")
                break
            adjust_learning_rate(model_optim, epoch+1, self.args)
        
        return self.model

    def test(self, cross,ii):
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
        print('test shape:', preds.shape, trues.shape)
        preds = preds.reshape(-1, preds.shape[-2], preds.shape[-1])
        trues = trues.reshape(-1, trues.shape[-2], trues.shape[-1])
        print('test shape:', preds.shape, trues.shape)

        folder_path = './results/' + self.args.setting +'/'
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        np.save(folder_path+'pred.npy', preds)
        np.save(folder_path+'true.npy', trues)

        mae, mse, rmse, mape, mspe, rse, corr = metric(preds, trues)
        print('mse:{}, mae:{},rmse{},mape{},mspe{},rse{},corr{}'.format(mse, mae, rmse, mape, mspe, rse, corr))
        np.save(folder_path+'metrics.npy', np.array([mae, mse, rmse, mape, mspe,rse,corr]))
       
        print("本次实验测试结果：",preds)
        print("真实值：",trues)

        return

    def predict(self, cross, ii, load=False):
        pred_data, pred_loader = self._get_data(cross,flag='pred')
        
        if load:
            path = os.path.join(self.args.checkpoints, self.args.setting)
            if self.args.data_fine_tune:
                best_model_path = path+'/'+self.args.data_fine_tune+str(cross)+'checkpoint.pth'
            else:
                best_model_path = path+'/'+'checkpoint.pth'
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
        print("本次实验预测未来结果：",preds)

        return

    def _process_one_batch(self, dataset_object, batch_x, batch_y, batch_label,batch_station):
        batch_x = batch_x.float().to(self.device)
        batch_y = batch_y.float()
        batch_label_onehot = F.one_hot(batch_label, num_classes=4).float().to(self.device)

        if self.args.padding==0:
            dec_inp = torch.zeros([batch_y.shape[0], self.args.pred_len, batch_y.shape[-1]]).float()
        elif self.args.padding==1:
            dec_inp = torch.ones([batch_y.shape[0], self.args.pred_len, batch_y.shape[-1]]).float()
        dec_inp = torch.cat([batch_y[:,:self.args.label_len,:], dec_inp], dim=1).float().to(self.device)

        if self.args.use_amp:
            with torch.cuda.amp.autocast():
                if self.args.output_attention:
                    outputs_predict, outputs_classify,attn = self.model(batch_x, dec_inp)[0,1,2]
                else:
                    outputs_predict, outputs_classify = self.model(batch_x, dec_inp)
        else:
            if self.args.output_attention:
                outputs_predict, outputs_classify,attn = self.model(batch_x, dec_inp)
            else:
                outputs_predict, outputs_classify = self.model(batch_x, dec_inp)
        if self.args.inverse:
            outputs_predict = dataset_object.inverse_transform(outputs_predict)
        f_dim = -1 if self.args.features=='MS' else 0
        batch_y = batch_y[:,-self.args.pred_len:,f_dim:].to(self.device)
        if self.args.output_attention:
            return outputs_predict, batch_y, outputs_classify, batch_label_onehot, attn
        else:
            return outputs_predict, batch_y, outputs_classify, batch_label_onehot

    def retest(self, cross, ii, flag):
        test_data, test_loader = self._get_data(cross, flag)
        torch.set_num_threads(4)

        path = os.path.join(self.args.checkpoints, self.args.setting)
        if self.args.retest_data_fine_tune:
            best_model_path = path+'/'+self.args.retest_data_fine_tune+str(cross)+'checkpoint.pth'
        else:
            best_model_path = path+'/'+str(cross)+'checkpoint.pth'
        state_dict = self._load_model(best_model_path)
        self.model.load_state_dict(state_dict)
        self.model.eval()
        
        preds = []
        trues = []
        pred_labels = []
        batch_labels = []
        attns = []
        batchs_x = []

        for i, (batch_x, batch_y, batch_label, batch_station) in enumerate(tqdm(test_loader)):

            if self.args.output_attention:
                pred, true, pred_label,true_label,attn = self._process_one_batch(
                    test_data, batch_x, batch_y, batch_label, batch_station)
            else:
                pred, true, pred_label,true_label = self._process_one_batch(
                    test_data, batch_x, batch_y, batch_label, batch_station)  

            preds.append(pred.detach().cpu().numpy())
            trues.append(true.detach().cpu().numpy())
            pred_labels.append(pred_label.detach().cpu().numpy())
            batch_labels.append(batch_label.detach().cpu().numpy())
            if self.args.output_attention:
                attns.append(attn[0].detach().cpu().numpy())
            batchs_x.append(batch_x.detach().cpu().numpy())

        folder_path = './results/' +self.args.setting +'/'
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        preds = np.array(preds)
        trues = np.array(trues)
        batchs_x = np.array(batchs_x)
        batchs_x = batchs_x.reshape(-1, batchs_x.shape[-2], batchs_x.shape[-1])
        preds = preds.reshape(-1, preds.shape[-2], preds.shape[-1])
        trues = trues.reshape(-1, trues.shape[-2], trues.shape[-1])

        pred_labels = np.array(pred_labels)
        batch_labels = np.array(batch_labels)
        pred_labels5 = pred_labels.reshape(-1, pred_labels.shape[-1])

        pred_labels = np.argmax(pred_labels5, axis=1)
        pred_labels = pred_labels.reshape(-1, 1)
        batch_labels = batch_labels.reshape(-1, 1)
        if self.args.output_attention:
            attns = np.array(attns)
        if self.args.output_attention:
            attns = attns.reshape(-1, attns.shape[-3], attns.shape[-2], attns.shape[-1])

       
        corr_sub = np.zeros((test_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        mse_sub = np.zeros((test_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        mae_sub = np.zeros((test_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        rmse_sub = np.zeros((test_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        rse_sub = np.zeros((test_data.data_num,preds.shape[1],preds.shape[2]),dtype = None, order = 'C')
        if self.args.output_attention:
            attn_sub = np.zeros((test_data.data_num,attns.shape[1],attns.shape[2],attns.shape[3]),dtype = None, order = 'C')

        acc_sub = np.zeros((test_data.data_num,pred_labels.shape[1]),dtype = None, order = 'C')
        pre_sub = np.zeros((test_data.data_num,pred_labels.shape[1]),dtype = None, order = 'C')
        rec_sub = np.zeros((test_data.data_num,pred_labels.shape[1]),dtype = None, order = 'C')
        f1_sub = np.zeros((test_data.data_num,pred_labels.shape[1]),dtype = None, order = 'C')
        label_sub = np.zeros((test_data.data_num,pred_labels.shape[1]),dtype = None, order = 'C')

        if self.args.data_path == "MULTI_RUMI_EMO.hdf5":
            labels=[0,1,2,3]
            labelss = ['rest', 'emotion', 'rumination', 'distraction']
        elif self.args.data_path == 'MDD_4task.hdf5':
            labels=[0,1,2,3]
            labelss = ["mdd","nc","bd","scz"]
        elif self.args.data_path == 'Rumination_Suzhou.hdf5':
            labels=[0,1,2,3]
            labelss = ['rest', 'emotion', 'rumination', 'distraction']

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

                acc_sub[sub] = ACC(pred_labels[corr_board0:corr_board1],batch_labels[corr_board0:corr_board1])
                pre_sub[sub] = PRE(pred_labels[corr_board0:corr_board1],batch_labels[corr_board0:corr_board1],labels)
                rec_sub[sub] = RECALL(pred_labels[corr_board0:corr_board1],batch_labels[corr_board0:corr_board1],labels)
                f1_sub[sub] = F1(pred_labels[corr_board0:corr_board1],batch_labels[corr_board0:corr_board1],labels)
                label_sub[sub] = batch_labels[corr_board0:corr_board1][0] 
            else:
                corr_sub[sub] = CORR(preds[corr_board0:],trues[corr_board0:])
                mse_sub[sub] = MSE(preds[corr_board0:],trues[corr_board0:])
                mae_sub[sub] = MAE(preds[corr_board0:],trues[corr_board0:])
                rmse_sub[sub] = RMSE(preds[corr_board0:],trues[corr_board0:])
                rse_sub[sub] = RSE(preds[corr_board0:],trues[corr_board0:])
                if self.args.output_attention:
                    attn_sub[sub] = attns[corr_board0:].mean(0)
                acc_sub[sub] = ACC(pred_labels[corr_board0:],batch_labels[corr_board0:])
                pre_sub[sub] = PRE(pred_labels[corr_board0:],batch_labels[corr_board0:],labels)
                rec_sub[sub] = RECALL(pred_labels[corr_board0:],batch_labels[corr_board0:],labels)
                f1_sub[sub] = F1(pred_labels[corr_board0:],batch_labels[corr_board0:],labels)
                label_sub[sub] = batch_labels[corr_board0:][0]

        corr_mean2 = corr_sub.mean(0).mean(0).mean(0)
        mse_mean2 = mse_sub.mean(0).mean(0).mean(0)
        mae_mean2 = mae_sub.mean(0).mean(0).mean(0)
        rmse_mean2 = mse_sub.mean(0).mean(0).mean(0)
        rse_mean2 = rse_sub.mean(0).mean(0).mean(0)

        acc_mean2 = acc_sub.mean(0)
        pre_mean2 = pre_sub.mean(0)
        rec_mean2 = rec_sub.mean(0)
        f1_mean2 = f1_sub.mean(0)
        corr_mean1 = CORR(preds,trues).mean(0).mean(0)
        mse_mean1 = MSE(preds,trues).mean(0).mean(0)
        mae_mean1 = MAE(preds,trues).mean(0).mean(0)
        rmse_mean1 = RMSE(preds,trues).mean(0).mean(0)
        rse_mean1 = RSE(preds,trues).mean(0).mean(0)

        acc_mean1 = ACC(pred_labels,batch_labels)
        pre_mean1 = PRE(pred_labels,batch_labels,labels)
        rec_mean1 = RECALL(pred_labels,batch_labels,labels)
        f1_mean1 = F1(pred_labels,batch_labels,labels)
        report = CLASS_REPORT(pred_labels,batch_labels,labels,labelss)
        confusion_matrix = CONFUSION_MATRIX(pred_labels,batch_labels)
        right_count = np.sum(acc_sub>0.5)/acc_sub.shape[0]
        LOGGER.info("corr_mean2:{}, mse_mean2:{}, mae_mean2:{}, rmse_mean2:{}, rse_mean2:{}".format(corr_mean2,mse_mean2,mae_mean2,rmse_mean2,rse_mean2))
        LOGGER.info("corr_mean1:{}, mse_mean1:{}, mae_mean1:{}, rmse_mean1:{}, rse_mean1:{}".format(corr_mean1,mse_mean1,mae_mean1,rmse_mean1,rse_mean1))
        LOGGER.info("acc_mean2:{}, pre_mean2:{}, rec_mean2:{}, f1_mean2:{}".format(acc_mean2,pre_mean2,rec_mean2,f1_mean2))
        LOGGER.info("acc_mean1:{}, pre_mean1:{}, rec_mean1:{}, f1_mean1:{}".format(acc_mean1,pre_mean1,rec_mean1,f1_mean1))
        LOGGER.info("report:{}".format(report))
        LOGGER.info("report:{}".format(confusion_matrix))
        LOGGER.info("right_count:{}".format(right_count))



        np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'repred.npy', preds) 
        np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'retrue.npy', trues) 
        np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'rebatch_labels.npy', batch_labels) 
        np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'repred_labels.npy', pred_labels) 
        np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'pred_labels_5.npy', pred_labels5) 

        np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'remetrics.npy', np.array([mse_sub, mae_sub, rmse_sub, rse_sub,corr_sub]))
        np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'remetrics_classify.npy', np.array([acc_sub, pre_sub, rec_sub, f1_sub]))
        np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'remetrics_mean.npy', np.array([mse_mean2, mae_mean2, rmse_mean2, rse_mean2, corr_mean2, acc_mean1, pre_mean1, rec_mean1, f1_mean1, acc_mean2, pre_mean2, rec_mean2, f1_mean2,right_count]))
        if self.args.output_attention:
            np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'attn.npy', attn_sub)
        np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'label_sub.npy', label_sub) 
        np.save(folder_path+self.args.data_path.replace('.hdf5', '')+'_'+flag+self.args.retest_data_fine_tune+str(cross)+'sub_info.npy', np.array([test_data.data_num, test_data.num_win])) 